# calls lme4 to estimate the parameters of a local distribution.
fitME = function(node, parents, group, data) {
  
  rhs = paste(c("1", parents), collapse = " + ")
  formula = paste(node, "~", rhs, "+", paste0("(", rhs, " | ", group, ")"))
  
  # Decide the family to use for the mixed model
  if (is.factor(data[[node]])) {
    model = try(glmer(formula, data = data, 
                      control = glmerControl(calc.derivs = FALSE,
                                             optCtrl=list(maxfun=3e5)),
                      family = binomial(link = "logit")
    ))
  } else {
    model = try(lmer(formula, data = data, 
                     control = lmerControl(calc.derivs = FALSE,
                                           optCtrl=list(maxfun=3e5))))
  }
  if (is(model, "try-error"))
    return(NULL)
  
  return(model)
  
}

scoreME =  function(node, parents, data, args = list()) {
  
  if (node == args$group)
    return(0)
  
  # mixed effects models with performance tuning to speed up learning.
  # require(lme4)
  # environment(nloptwrap)$defaultControl$maxeval = 5 * 1e3
  # environment(nloptwrap)$defaultControl$xtol_abs = 1e-6
  # environment(nloptwrap)$defaultControl$ftol_abs = 1e-6
  
  # result of fitME is lme, so we can take the BIC of that
  local.distribution = fitME(node = node, 
                             parents = setdiff(parents, args$group),
                             group = args$group, data = data)
  
  if (is.null(local.distribution))
    return(-Inf)
  
  return(-BIC(local.distribution))
  
}

# calls lme4 to estimate the parameters of a local distribution.
fitPooled = function(node, parents, data) {
  
  rhs = paste(c("1", parents), collapse = " + ")
  formula = paste(node, "~", rhs)
  
  print(formula)
  # Decide the family to use for the mixed model
  if (is.factor(data[[node]])) {
    model = try(glm(formula, data = data, family = binomial(link = "logit")))
  } else {
    model = try(lm(formula, data = data))
  }
  
  if (is(model, "try-error"))
    return(NULL)
  
  return(model)
  
}

scorePooled =  function(node, parents, data, args = list()) {
  
  
  # mixed effects models with performance tuning to speed up learning.
  # require(lme4)
  # environment(nloptwrap)$defaultControl$maxeval = 5 * 1e3
  # environment(nloptwrap)$defaultControl$xtol_abs = 1e-6
  # environment(nloptwrap)$defaultControl$ftol_abs = 1e-6
  
  # result of fitME is lme, so we can take the BIC of that
  local.distribution = fitPooled(node = node, parents = parents, data = data)
  
  if (is.null(local.distribution))
    return(-Inf)
  
  return(-BIC(local.distribution))
  
}

 


ldistME = function(node, parents, group, data) {
  
  # mixed effects models with performance tuning to speed up learning.
  # environment(nloptwrap)$defaultControl$maxeval = 5 * 1e3
  # environment(nloptwrap)$defaultControl$xtol_abs = 1e-6
  # environment(nloptwrap)$defaultControl$ftol_abs = 1e-6
  
  model = fitME(node = node, parents = parents, group = group, data = data)
  ngroups = data[[group]] %>% unique %>% length
  
  ldist = list(
    coef = array(0, dim = c(length(parents) + 1, ngroups),
                 dimnames = list(c("(Intercept)", parents), NULL)),
    sd = rep(0, ngroups))
  
  if (is.null(model)) {
    
    # retry with more iterations and a more forgiving precision.
    # environment(nloptwrap)$defaultControl$maxeval = 1e5
    # environment(nloptwrap)$defaultControl$xtol_abs = 1e-5
    # environment(nloptwrap)$defaultControl$ftol_abs = 1e-5
      
    model = fitME(node = node, parents = parents, group = group, data = data)
      
    # reset performance tuning.
    # environment(nloptwrap)$defaultControl$maxeval = 5 * 1e3
    # environment(nloptwrap)$defaultControl$xtol_abs = 1e-6
    # environment(nloptwrap)$defaultControl$ftol_abs = 1e-6
    
    if (is.null(model)) {
      
      warning("failed to fit node ", node, " with lme4.")
      return(ldist)
    }
  
  }
  
  # add the fixed effects to all columns.
  fixefs = fixef(model)
  for (i in seq(ncol(ldist$coef)))
    ldist$coef[, i] = fixefs
  
  # update each column with the random effects from the corresponding group.
  ranefs = ranef(model)[[group]]
  for (i in seq(ncol(ldist$coef)))
    ldist$coef[, i] = ldist$coef[, i] + as.numeric(ranefs[i, ])
  
  # compute the standard errors for the residuals in the groups.
  resids = resid(model)
  g = data[[group]] %>% unique
  for (i in seq(ngroups))
    ldist$sd[i] = sd(resids[data[[group]] == g[i]])
  
  # replace NA coefficients and standard errors with zeroes to avoid errors in custom.fit().
  ldist$coef[is.na(ldist$coef)] = 0
  ldist$sd[is.na(ldist$sd)] = 0
  
  return(ldist)
                
}

ldistLM = function(node, data) {
  
  formula = paste(node, "~", 1)
  model = lm(formula, data = data)
  
  list(
    coef = model$coef,
    sd = sd(model$residuals)
  )
}

rlmebnfit = function(nnodes, ngroups, prob, homogeneous = FALSE) {
  
  nodes = paste0("X", seq(nnodes))
  graph = random.graph(nodes, num = 1, method = "ordered", prob = prob)
  
  while(igraph::is_connected(as.igraph(graph)) == FALSE)
    graph = random.graph(nodes, num = 1, method = "ordered", prob = prob)
  
  local.distributions = structure(vector(nnodes, mode = "list"), names = nodes)
  
  for (node in nodes) {
    par = parents(graph, node)
    betas = matrix(rep(0, ngroups * (length(par) + 1)),
                   nrow = length(par) + 1, ncol = ngroups,
                   dimnames = list(c("(Intercept)", par), seq(ngroups) - 1))
    
    if (homogeneous) {
      for (p in c("(Intercept)", par))
        betas[p, ] = rnorm(1, mean = rnorm(1, mean = 2, sd = 1),
                           sd = sqrt(rchisq(1, 1)))
    } else {
      for (p in c("(Intercept)", par))
        betas[p, ] = rnorm(ngroups, mean = rnorm(1, mean = 2, sd = 1),
                           sd = sqrt(rchisq(1, 1)))
    }
    
    sds = rep(sqrt(rchisq(1, 1)), ngroups)
    
    local.distributions[[node]] = list(coef = betas, sd = sds)
    
  }
  
  graph = add.node(graph, "G")
  for (node in nodes)
    graph = set.arc(graph, from = "G", to = node)
  
  local.distributions$G = matrix(rep(1/ngroups, ngroups), 
                                 ncol = ngroups,
                                 dimnames = list(NULL, paste0("g", seq(ngroups))))
  
  
  bn = custom.fit(graph, local.distributions)
  
  for (node in setdiff(node.ordering(bn), "G")) {
    
    pp = setdiff(parents(bn, node), "G")
    
    if (length(pp) == 0)
      next
    else
      ff = paste(node, "~", paste(pp, collapse = "+"))
    
    data = rbn(bn, 100 * nparams(bn))
    
    for (g in levels(data$G)) {
      model = lm(as.formula(ff), data[data$G == g, ])
      exp.var = var(resid(model)) / var(data[data$G == g, node])
      local.distributions[[node]]$sd[levels(data$G) == g] = sqrt(var(fitted(model)) / 0.85 * 0.15)
      
    }
    
    bn = custom.fit(graph, local.distributions)
    
  }
  
  return(bn)
  
}