##############################################################
# Kernel function with its respective growth*survival kernels
# Will Vieira
# January 20, 2019
##############################################################



######################
 # 1. Growth
 # 2. Survival (1 - mortality)
 # 3. Growth * survival kernel
 # 3. Full kernel
 ######################



# Variable growth function
growth_Var = function(y, x, Var, ...)
{
  mean = growth_XEC(x, E, C, parsGrowth)
  growth = dgamma(y, shape = mean^2/Var, rate = mean/Var)
  growth[is.infinite(growth)] = 0 # TODO: clean this dirty code

  return( growth )
}

# Variable survival function derived from the deterministic mortality function
survival = function(x, ...)
{
  surv = 1 - mort_XEC(x, E, C, parsMort)

  return( surv )
}



## Survival/growth kernel
P_xEC = function(y, x, Var, E, C, parsMort, parsGrowth)
{
  pkernel = growth_Var(y-x, x, Var, E, C, parsGrowth) *
            survival(x, E, C, parsMort)

  return( pkernel )
}



mkKernel = function(m, U, L, E, C, parsMort, parsGrowth, parsFec)
{
  # mesh points
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h

  P <- h * outer(meshpts, meshpts, P_xEC, E, C, parsMort, parsGrowth, Var = 3)
  F <- h * outer(meshpts, meshpts, f.yx, parsFec)

  K <- P + F

  return(list(K = K, meshpts = meshpts, P = P, F = F))
}
