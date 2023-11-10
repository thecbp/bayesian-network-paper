source("helpers.R")

# Generate the parameter permutations
# Edited on my personal computer to speed this up
params_v2 = expand.grid(
  n.subj = c(20),
  N = c(30),
  avg.parents = c(2),
  sampsizes = c(60)
) %>% as_tibble() %>% 
  mutate(
    file = pmap(list(n.subj, N, avg.parents, sampsizes), function(ns, n, ap, ss) {
      paste0("F", ns, "-N", n, "-A", ap, "-samp-", ss)
    }) %>% unlist
  ) %>% 
  arrange(N, n.subj, avg.parents, sampsizes) %>% 
  mutate(
    n = c(0)
  )

currently_done = list.files("/Volumes/T7/bayesian network-paper/simluations-v2") %>% 
  str_replace("-iter-[:digit:].rds", "") %>% 
  str_replace("metrics-", "") %>% 
  unique

# Check to see how many runs have been done for each
check = list.files("/Volumes/T7/bayesian network-paper/simluations-v2") %>% 
  str_replace("-iter-[:digit:].rds", "") 
done = tibble(
  perm = check
) %>% 
  group_by(perm) %>% summarize(n = n()) %>% 
  mutate(
    file = str_replace(perm, "metrics-", "") %>% unlist()
  )

params_v2 = params_v2 %>% 
  left_join(done, by = "file") %>% 
  replace_na(list(n = 0))
  

for (row in 1:nrow(params_v2)) {
  
  if (params_v2$n[row] == 5) {
    done_log = paste0(lubridate::now(), params_v2$file[row], " is done, skipping")
    write(done_log, file = "simulation-log-v2.txt", append = TRUE, sep = "\n")
    next
  }
  
  seed = 1 
  successes = params_v2$n[row]
  run = params_v2$n[row]
  
  param_log = paste0(lubridate::now(), ": Starting ", paste0("F", params_v2$n.subj[row],
                                                             "-N", params_v2$N[row],
                                                             "-A", params_v2$avg.parents[row],
                                                             "-samp-", params_v2$sampsizes[row]))
  write(param_log, file = "simulation-log-v2.txt", append = TRUE, sep = "\n")
  
  while (successes < 5) {
    
    run = run + 1
    
    run_log = paste0(lubridate::now(), ": Starting run ", run, " for ", 
                       paste0("F", params_v2$n.subj[row],
                              "-N", params_v2$N[row],
                              "-A", params_v2$avg.parents[row],
                              "-samp-", params_v2$sampsizes[row]))
    write(run_log, file = "simulation-log-v2.txt", append = TRUE, sep = "\n")
    
    
    sim = generate_abn_data(n.subj = params_v2$n.subj[row], 
                            N = params_v2$N[row], 
                            avg.parents = params_v2$avg.parents[row],
                            samp.size = params_v2$sampsizes[row],
                            prop.gauss = 0.5,
                            prop.pos.fixef = 0.5,
                            seed = seed)
    
    
    print(paste0("Parameters: ", paste0("metrics-F", params_v2$n.subj[row],
                                        "-N", params_v2$N[row],
                                        "-A", params_v2$avg.parents[row],
                                        "-samp-", params_v2$sampsizes[row],
                                        ". Run: ", run )))
    
    test = tryCatch(
      expr = {
        simulation_suite(sim)
      }, 
      error = function(e){         
        message("Bad seed.")
      })
    
    # Save the file if it succeeds
    if (!is.null(test)) {
      
      successes = successes + 1 
      
      file = paste0("metrics-F", params_v2$n.subj[row],
                    "-N", params_v2$N[row],
                    "-A", params_v2$avg.parents[row],
                    "-samp-", params_v2$sampsizes[row], 
                    "-iter-", successes, 
                    ".rds")
      
      saveRDS(test, file)
      
      seed = seed + 1
      
      
    } else {
      
      fail_log = paste0(lubridate::now(), ": Failed run at #", run, ". New seed.")
      write(fail_log, file = "simulation-log-v2.txt", append = TRUE, sep = "\n")
      
      # Otherwise, increment the seed
      seed = seed + 1
    }
    
    
    
  }
  
  param_log = paste0(lubridate::now(), ": Ending ", paste0("F", params_v2$n.subj[row],
                                                             "-N", params_v2$N[row],
                                                             "-A", params_v2$avg.parents[row],
                                                             "-samp-", params_v2$sampsizes[row]))
  write(param_log, file = "simulation-log-v2.txt", append = TRUE, sep = "\n")
  
  
}






#########################

# DEBUGGING AREA / ERROR CATCHING

#########################

# samp 120, seed 5 causes the error in backward selection
# 01-05-2023: F10-N50-A1-samp-30: 
  # debugging report: some models accidentally use "1" as the reference, causing missnaming of parameters


sampsize = 30
seed = 1

file = paste0("metrics-F", cur$n.subj,
              "-N", cur$N,
              "-A", cur$avg.parents,
              "-samp-", sampsize, 
              "-iter-", 5, 
              ".rds")
print(file)

sim = generate_abn_data(n.subj = 10, 
                        N = 50, 
                        avg.parents = 1,
                        samp.size = sampsize,
                        seed = seed,
                        prop.gauss = 0.5,
                        prop.pos.fixef = 0.5)

 # Construct a naive BN
time = lubridate::now()
naive.bn = hc(sim$data, maxp = 5)
print(paste0("ME network started at: ", time, " and took ", lubridate::now() - time))

# Construct a mixed BN
time = lubridate::now()
me.bn = hc(sim$data, score = "custom", fun = scoreME, args = list(group = "G"), maxp = 5)
print(paste0("ME network started at: ", time, " and took ", lubridate::now() - time))

# Do stepwise selection for all variables
time = lubridate::now()
bw = build_stepwise_models(sim$data, sim$dists, "backward")
print(paste0("Backward selection started at: ", time, " and took ", lubridate::now() - time))

time = lubridate::now()
fw = build_stepwise_models(sim$data, sim$dists, "forward")
print(paste0("Forward selection started at: ", time, " and took ", lubridate::now() - time))

naive.dag = bn2adjmat(naive.bn)
me.dag = bn2adjmat(me.bn)

naive.missing.edges = sum((sim$dag - naive.dag) == 1)
me.missing.edges = sum((sim$dag - me.dag) == 1)
bw.missing.edges = sum((sim$dag - bw$adjmat) == 1)
fw.missing.edges = sum((sim$dag - fw$adjmat) == 1)

# extra edges
naive.extra.edges = sum((sim$dag - naive.dag) == -1)
me.extra.edges = sum((sim$dag - me.dag) == -1)
bw.extra.edges = sum((sim$dag - bw$adjmat) == -1)
fw.extra.edges = sum((sim$dag - fw$adjmat) == -1)

# correct edges, but flipped
naive.matchflip.edges = sum(sim$dag * t(naive.dag))
me.matchflip.edges = sum(sim$dag * t(me.dag))
bw.matchflip.edges = sum(sim$dag * t(bw$adjmat))
fw.matchflip.edges = sum(sim$dag * t(fw$adjmat))

# correct edge + orientation correct
naive.matched.edges = sum(sim$dag * naive.dag)
me.matched.edges = sum(sim$dag * me.dag)
bw.matched.edges = sum(sim$dag * bw$adjmat)
fw.matched.edges = sum(sim$dag * fw$adjmat)

# hamming distance: 
naive.hamming = naive.missing.edges + naive.extra.edges + naive.matchflip.edges 
me.hamming = me.missing.edges + me.extra.edges + me.matchflip.edges 
bw.hamming = bw.missing.edges + bw.extra.edges + bw.matchflip.edges 
fw.hamming = fw.missing.edges + fw.extra.edges + fw.matchflip.edges 



out = list()
out[["simulation"]] = sim
out[["true.dag"]] = sim$dag
out[["naive_bn"]] = naive.bn
out[["me_bn"]] = me.bn
out[["backward"]] = bw
out[["forward"]] = fw

out[["naive.missing.edges"]] = naive.missing.edges
out[["me.missing.edges"]] = me.missing.edges
out[["bw.missing.edges"]] = bw.missing.edges
out[["fw.missing.edges"]] = fw.missing.edges

out[["naive.extra.edges"]] = naive.extra.edges
out[["me.extra.edges"]] = me.extra.edges
out[["bw.extra.edges"]] = bw.extra.edges
out[["fw.extra.edges"]] = fw.extra.edges

out[["naive.matchflip.edges"]] = naive.matchflip.edges
out[["me.matchflip.edges"]] = me.matchflip.edges
out[["bw.matchflip.edges"]] = bw.matchflip.edges
out[["fw.matchflip.edges"]] = fw.matchflip.edges

out[["naive.matched.edges"]] = naive.matched.edges
out[["me.matched.edges"]] = me.matched.edges
out[["bw.matched.edges"]] = bw.matched.edges
out[["fw.matched.edges"]] = fw.matched.edges

out[["naive.hamming"]] = naive.hamming
out[["me.hamming"]] = me.hamming
out[["bw.hamming"]] = bw.hamming
out[["fw.hamming"]] = fw.hamming

saveRDS(out, file)
