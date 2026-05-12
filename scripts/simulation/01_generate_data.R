#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1 Model parameters ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(geomix)

## 1.1 K=8 parameter set (full simulation study) ----
params <- list(K = 8, alpha = matrix(c(10,17,14,19,7,11,6,12,
                      0.2,0.05,0.2,0.05,0.1,0.15,0.2,0.15),ncol=2),
     sigma2 = c(5,4,6,5,10,4,3,7),
     tau = 0.5,
     lL = c(0.5,0.5,0.6,0.4,0.6,0.5,0.7,0.3),
     lD = c(1,0.8,0.6,0.7,0.6,1,0.7,0.5)+1,
     beta = 0.6, h = 2,
     deltaK = 0.1, gammaScale = 0.5, kappa = 0.99)

## 1.2 K=3 parameter set (small demo) ----
params3 <- list(K = 3, alpha = matrix(c(10,19,7,0.1,0.15,0.2),ncol=2),
               sigma2 = c(5,4,7),
               tau = 1,
               lL = c(1,0.75,1.5),
               lD = c(2,3,4),
               beta = 0.6, h = 2,
               deltaK = 0.1, gammaScale = 0.5, kappa = 0.95)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2 Synthetic data generator ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 2.1 generate_geomix_synthetic: simulate a GeoMix dataset ----
#
# Arguments:
#   params    - list of model parameters (K, alpha, sigma2, tau, lL, lD, beta, h, ...)
#   nd/nx/ny  - depth, x, y grid dimensions
#   prop_obs  - proportion of lateral locations with Z2 observations
#   seed      - RNG seed
generate_geomix_synthetic <- function(params,nd = 20, nx = 20, ny = 20,
                                      prop_obs = 0.4, seed = 2025) {
  set.seed(seed)
  require(mvtnorm)
  list2env(params, envir = environment())

  ## 2.2 Helper: enforce monotone stratigraphy along depth ----
  enforce_monotone <- function(x) {
    n <- length(x)
    out <- x
    change_pts <- which(diff(x) < 0)
    bad1 <- 1 %in% change_pts
    if(bad1) out[1] <- out[2]
    change_pts <- change_pts[change_pts != 1]
    m <- length(change_pts)
    while(m > 0){
      check <- sum(sign(diff(out[change_pts[1] + c(-1,0,1)])))+1
      if(!(check %in% c(0,1))){
        cat(paste0("CHECK x: ",paste0(x,collapse = " "),"\n"))
      }
      if(check) out[change_pts[1]] <- out[change_pts[1]-1] else out[change_pts[1]+1] <- out[change_pts[1]]
      change_pts <- which(diff(out) < 0)
      bad1 <- 1 %in% change_pts
      if(bad1) out[1] <- out[2]
      change_pts <- change_pts[change_pts != 1]
      m <- length(change_pts)
    }
    return(out)
  }

  ## 2.3 Grid and neighbourhood structure ----
  coords <- expand.grid(d = 1:nd, x = 1:nx, y = 1:ny)
  coords$ID <- 1:nrow(coords)
  N1 <- nrow(coords)
  D <- nd
  L <- nx * ny
  neighbours <- getNeighbours(dims = c(nd,nx,ny), diagonals = 0)
  neighbourNum <- apply(neighbours>0,1,sum)
  weights <- matrix(1, nrow = nrow(coords), ncol = ncol(neighbours))
  penalty <- as.matrix(dist(1:K))

  ## 2.4 Simulate latent stratigraphy Y1 via Potts model ----
  layer_bounds <- floor(seq(1, nd + 1, length.out = K + 1))
  Y1 <- rep(0, N1)
  for (k in 1:K) {
    d_inds <- which(coords$d >= layer_bounds[k] & coords$d < layer_bounds[k + 1])
    Y1[d_inds] <- k
  }
  for(i in 1:10) Y1 <- samplePotts(beta=beta,Y1,K=K,weights,neighbours,neighbourNum,penalty)
  Y1 <- apply(array(Y1,dim = c(nd,nx,ny)),2:3,enforce_monotone)
  Y1[1,,] <- 1
  Y1 <- c(Y1)

  ## 2.5 Depth and lateral indices ----
  coords$locID <- with(coords, (y - 1) * nx + x)
  coords$dID <- coords$d
  locID <- coords$locID
  dID <- coords$dID

  loc_coords <- unique(coords[, c("x", "y")])
  distL <- as.matrix(dist(loc_coords))^2
  distD <- as.matrix(dist(1:D))^2

  depth <- coords$d
  X <- model.matrix(~1+depth)

  ## 2.6 Simulate latent geotechnical property Y2 ----
  Y2 <- numeric(N1)
  mu <- numeric(N1)
  sig <- numeric(N1)
  for (k in 1:K) {
    inds <- which(Y1 == k)
    if (length(inds) == 0) next
    Xk <- X[inds, , drop = FALSE]
    mk <- Xk %*% alpha[k, ]
    r2 <- distL[locID[inds], locID[inds]] / lL[k] +
      distD[dID[inds], dID[inds]] / lD[k]
    covk <- sigma2[k] * (1 + sqrt(3)*sqrt(r2)) * exp(-sqrt(3)*sqrt(r2))
    Y2[inds] <- as.numeric(rmvnorm(1, mk, sigma = (covk)))
    mu[inds] <- mk
    sig[inds] <- sqrt(sigma2[k])
  }

  ## 2.7 Subsample Z2 observations ----
  loc_select <- sort(sample(unique(locID), size = floor(prop_obs * length(unique(locID)))))
  Z2_ind <- which(locID %in% loc_select)
  N2 <- length(Z2_ind)
  Z2_nind <- setdiff(1:N1,Z2_ind)
  Z2 <- Y2[Z2_ind] + rnorm(N2,mean=0,sd=tau)
  Z2_test <- Y2[Z2_nind]

  ## 2.8 Simulate ground model stratigraphy Z1 from confusion model ----
  gamma0 <- exp(-as.matrix(dist(1:K)))+deltaK
  diag(gamma0) <- 0
  gamma <- gammaScale*t(apply(gamma0, 1, function(x) x / sum(x)))
  diag(gamma) <- 0

  locStack <- coords %>%
    group_by(locID) %>%
    summarise(minD = min(ID),
              maxD = max(ID),.groups="drop")
  Z1 <- numeric(N1)
  phimat <- pnorm(-0.5 * (sqrt(distD[-1, ]) + sqrt(distD[-nd, ])) / h)

  for(i in 1:nrow(locStack)){
    pis <- computePiStack(loc0 = locStack$minD[i], loc1 = locStack$maxD[i],
                          Y1Stack = Y1[locStack$minD[i]:locStack$maxD[i]],
                          gammaMat = gamma, phiStack = phimat, K= K,
                          kappa  = kappa)
    Z1[locStack$minD[i]:locStack$maxD[i]] <-
      sapply(1:nrow(pis),function(j) sample(1:K,1,prob=pis[j,]))
  }
  Z1 <- c(apply(array(Z1,dim = c(nd,nx,ny)),2:3,enforce_monotone))
  table(Y1,Z1)

  coords$group <- coords$d

  data <- data.frame(coords, Y1, Y2, mu, sig, Z1, Z2 = NA)
  data$Z2[Z2_ind] <- Z2

  return(list(
    K = K,
    data = data,
    dims = c(nd,nx,ny),
    lattice = cbind(coords,Z1),
    X = X,
    alpha = alpha,
    sigma2 = sigma2,
    tau = tau,
    lL = lL,
    lD = lD,
    beta = beta,
    diagonals = 0,
    penalty = penalty,
    gamma = gamma,
    h = h,
    Z2_test = Z2_test,
    test_ind = Z2_nind
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3 Generate and save ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir.create("data/simulation", recursive = TRUE, showWarnings = FALSE)

## 3.1 K=8 dataset for simulation study ----
simulatedK8 <- generate_geomix_synthetic(params = params, prop_obs = 0.4)
saveRDS(simulatedK8, 'data/simulation/synthetic_data.rds')
message("Saved data/simulation/synthetic_data.rds")

## 3.2 K=3 dataset (small demo / misspecification figure) ----
simulatedK3 <- generate_geomix_synthetic(params = params3, prop_obs = 0.2)
saveRDS(select(simulatedK3$data, ID, d, x, y, locID, Z1, Z2), 'data/simulation/offshore.rds')
message("Saved data/simulation/offshore.rds")
