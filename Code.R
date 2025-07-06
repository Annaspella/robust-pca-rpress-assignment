n <- nrow(data)
p <- ncol(data)

scale=FALSE
signflip=TRUE 
scale <- FALSE
signflip <- TRUE 
## ___Step 1___: Reduce the data space to the affine subspace spanned by the n observations
##  Apply svd() to the mean-centered data matrix. If n > p we use the kernel approach -
##  the decomposition is based on computing the eigenvalues and eigenvectors of(X-m)(X-m)'
# Xsvd <- kernelEVD(data, scale=scale, signflip=signflip)
Xsvd <- classPC(data, scale=scale, signflip=signflip)

if (Xsvd$rank == 0)
  stop("All data points collapse!")

## VT::27.08.2010: introduce 'scale' parameter; return the scale in the value object
##
myscale <- vector('numeric', p) + 1
if (scale)
  myscale <- Xsvd$scale

##
## verify and set the input parameters: alpha, k and kmax
## determine h based on alpha and kmax, using the function h.alpha.n()
##

## VT::06.11.2012 - kmax <= floor(n/2) is too restrictive
##    kmax <- max(min(floor(kmax), floor(n/2), Xsvd$rank),1)
kmax <- max(min(floor(kmax), Xsvd$rank), 1)


if ((k <- floor(k)) < 0)
  k <- 0
else if (k > kmax) {
  warning(paste("The number of principal components k = ", k, " is larger than kmax = ", kmax, "; k is set to ", kmax,".", sep=""))
  k <- kmax
}

if (missing(alpha))
{
  default.alpha <- alpha
  if (is.null(h)) {
    h <- min(ceiling(alpha*n)+1, n)
  }
  alpha <- h/n
  if (h > n) {
    alpha <- default.alpha
    h <- ceiling(alpha*n)+1
    warning(paste("h should be smaller than n = ", n, ". It is set to its default value ", h, ".", sep=""))
  }
} else
{
  
  if (alpha < 0.5 | alpha > 1)
    stop("Alpha is out of range: should be between 1/2 and 1")
  
  hdef <- h
  h <- ceiling(alpha*n)+1
  
  if (!is.null(hdef))
    warning(paste("Both alpha and h are given, h is set to its default value ", h, ".", sep=""))
}

T1 <- scale(data, center=Xsvd$center, scale=Xsvd$scale) %*% Xsvd$loadings
center <- Xsvd$center
rot <- Xsvd$loadings
##
## ___Step 2___: Always apply full ROBPCA (no MCD if p<<n)
##

if (trace)
  cat("\nApplying the projection method of Hubert.\n")

##
##  For each direction 'v' through 2 data points we project the n data points xi on v
##  and compute their robustly standardized absolute residual
##  |xi'v - tmcd(xj'v)|/smcd(xj'v)
##

if (ndir=="all" & n>500) {
  warning("Computing all directions for this value of n can take very long, we
            recommend to set ndir to 5000.")
}

if (skew) {
  # Adjusted outlyingness from mrfDepth package, fast enough to compute all directions
  outl_obj <- adjOutl(T1, options = list(type="Rotation", ndir=ndir))
} else {
  # Outlyingness from mrfDepth package, fast enough to compute all directions
  outl_obj <- outlyingness(T1, options = list(type="Rotation", h=h, ndir=ndir, stand="unimcd"))
}

outl <- outl_obj$outlyingnessX
if (!is.numeric(outl)) {
  stop("Something went wrong with the outlyingness computations.")
}

if (any(!is.finite(outl)) | any(is.nan(outl))) {
  stop("The obtained outlyingness is not finite or NaN.")
}

if (min(outl)<10^(-16)) {
  stop("The obtained outlyingness is not strictly positive.")
}

if (length(outl)!=nrow(T1)) {
  stop("The obtained outlyingness does not have the correct length.")
}

@@ -201,8 +201,7 @@
  
  kmax <- min(Xh.svd$rank, kmax)
if (k == 0) {
  test <- which(Xh.svd$eigenvalues/Xh.svd$eigenvalues[1] <= 
                  0.001)
  test <- which(Xh.svd$eigenvalues/Xh.svd$eigenvalues[1] <= 0.001)
  k <- if (length(test) != 0) 
    min(min(Xh.svd$rank, test[1]), kmax)
  else min(Xh.svd$rank, kmax)
  cumulative <- cumsum(Xh.svd$eigenvalues[1:k])/sum(Xh.svd$eigenvalues)
  if (cumulative[k] > 0.8) {
    k <- which(cumulative >= 0.8)[1]
  }
  if (trace) 
    cat(paste("The number of principal components is set by the algorithm. It is set to ", 
              k, ".\n", sep = ""))
}
if (trace)
  cat("\nEigenvalues: ", Xh.svd$eigenvalues, "\n")

## perform extra reweighting step
if (k != Xsvd$rank)
{
  ## VT::27.08.2010 - bug report from Stephen Milborrow: if n is small relative to p
  ## k can be < Xsvd$rank but larger than Xh.svd$rank - the number of observations in
  ## Xh.svd is roughly half of these in Xsvd
  k <- min(Xh.svd$rank, k)
  
  XRc <- sweep(T1, 2, Xh.svd$center,'-')
  Xtilde <- XRc %*% Xh.svd$loadings[,1:k] %*% t(Xh.svd$loadings[,1:k])
  Rdiff <- XRc - Xtilde
  odh <- apply(Rdiff, 1, vecnorm)
  
  #Cutoff for the ODs
  tmp <- coOD(od=odh, h=h, skew=skew)
  odh <- tmp$od
  cutoffodh <- tmp$co.od
  
  H1 <- (odh <= cutoffodh)
  #Xh.svd <- classSVD(T1[H1,])
  Xh.svd <- classPC(T1[H1,])
  k <- min(Xh.svd$rank, k)
} else {
  H1 <- H0
}
return(list(H0=H0, H1=H1, k=k, alpha=alpha, h=h))
}
rospca_part2 <- function(X, H0, H1, k, kmax=10, alpha=0.75, h=NULL, grid=TRUE, lambda=10^(-6), stand=TRUE, 
                         sparse="penalty", para, skew=FALSE, ndir="all") {
  
  d <- dim(X); n <- d[1]; p <- d[2]
  X <- as.matrix(X) # Make matrix out of data matrix
  if (is.null(h)) h <- floor(alpha*n)
  
  
  # Obtain the observations from H1
  Xh1 <- X[H1,]
  
  
  ############################
  ## Sparse part
  ############################
  
  Xh1_stand <- Xh1
  X_H1stand <- X
  
  mu <- apply(X[H1,], 2, median)
  
  if (stand) {
    scale <- apply(X[H1,], 2, qn)
  } else {
    scale <- rep(1, p)
  }
  
  
  # Standardise the data 
  X_H1stand <- matStandVec(X, mu, scale)
  Xh1_stand <- matStandVec(Xh1_stand, mu, scale)
  D <- diag(scale)
  
  
  
  if (grid) {
    # SCoTLASS on the standardised regular observations
    P1 <- unclass(SPcaGrid(Xh1_stand, k=k, lambda=lambda, method="sd", scale=FALSE, center=rep(0,p),
                           glo.scatter=1, maxiter=75)@loadings)
  } else {
    # SPCA on the standardised regular observations
    P1 <- spca(Xh1_stand, K=k, para=para, sparse=sparse, lambda=lambda)$loadings
  }
  
  # Extra sparse reweighting step
  # Look for non-zero rows in the loadings matrix
  # index=which(rowSums(abs(P1))>10^(-10))
  # 1 if at least 1 non-zero element
  f <- function(x) {
    any(abs(x)>10^(-10))
  }
  index <- which(apply(P1, 1, f))
  Xtilde <- X_H1stand%*%P1%*%t(P1)
  # Remove variables that have zero loadings on all PCs, when computing ODs
  odh <- apply(X_H1stand[,index] - Xtilde[,index], 1, vecnorm)
  
  tmp <- coOD(od=odh, h=h, skew=skew)
  odh <- tmp$od
  cutoffodh <- tmp$co.od
  
  H2 <- (odh <= cutoffodh)
  
  Xh2_stand <- X_H1stand[H2,]
  
  P2 <- matrix(0,p,k)
  if (grid) {
    pp <- length(index)
    # SCoTLASS on the standardised observations with indices in H2
    P2[index,] <- SPcaGrid(Xh2_stand[,index], k=k, lambda=lambda, method="sd", scale=FALSE, center=rep(0,pp),
                           glo.scatter=1, maxiter=75)@loadings
  } else {
    # SPCA on the standardised observations with indices in H2
    P2[index,] <- spca(Xh2_stand[,index], K=k, para=para, sparse=sparse, lambda=lambda)$loadings
  }
  
  TT <- sweep(X, 2, mu, "-")%*%P2 # Scores matrix (now we use X because scores of all observations!)
  
  
  # Eigenvalues are variances of sparse PCs of observations in H2
  # We use a robust measure of scale because good leverage points can still be present.
  l1 <- apply(TT[H2,], 2, qn)^2   
  
  ##
  # SD reweighting
  
  
  
  if (skew) {
    
    # Score distances using adjusted outlyingness
    outl2 <- adjOutl(TT, options = list(type="Rotation", ndir=ndir))$outlyingnessX
    
    # Set H3 of all observations in H2 with sd <= cutoff
    H3 <- rep(FALSE, n)
    H3[H2] <- (outl2[H2] < coOD(outl2[H2], skew=TRUE)$co.od)
    H3 <- as.logical(H3)
    
  } else {
    
    # Compute score distances of observations in H2
    A <- distPCA(X=X[H2,], Tn=TT[H2,], P=P2, l=l1, mu=mu, h=sum(H2), skew=FALSE, co.od=FALSE)
    
    # Set H3 of all observations in H2 with sd <= cutoff
    H3 <- rep(FALSE, n)
    H3[H2] <- (A$sd<A$cutoff.sd)
    H3 <- as.logical(H3)
  }
  
  
  # New centre
  mu <- colMeans(X[H3,])
  
  # New scores
  TT <- sweep(X, 2, mu, "-")%*%P2 #Scores matrix (now we use X because scores of all observations!)
  
  # New eigenvalues
  l1 <- apply(TT[H3,], 2, var)   
  # End SD reweighting
  
  # Sort eigenvalues in decreasing order
  sorted <- sort(l1, index=TRUE, decreasing=TRUE)
  l <- sorted$x; ix <- sorted$ix
  
  P <- P2[,ix] # Final loadings matrix
  Tn <- TT[,ix] # Final scores matrix
  
  # Change names of columns
  s <- paste("PC", 1:k, sep="")
  
  names(l) <- s
  colnames(P) <- s 
  colnames(Tn) <- s 
  colnames(D) <- colnames(X)
  # Change names of rows
  rownames(P) <- colnames(X)
  rownames(D) <- colnames(X)
  if (!is.null(rownames(X))) {
    rownames(Tn) <- rownames(X)
  } else {
    rownames(Tn) <-  1:n
  }
  names(H0) <- rownames(Tn)
  names(H1) <- rownames(Tn)
  names(H2) <- rownames(Tn)
  names(H3) <- rownames(Tn)
  names(mu) <- colnames(X)
  
  
  # Compute distances and cutoffs
  A <- distPCA(X=X, Tn=Tn, P=P, l=l, mu=mu, h=h, skew=skew)
  sd <- A$sd
  od <- A$od
  cutoff.sd <- A$cutoff.sd
  cutoff.od <- A$cutoff.od
  
  # Adjust score distances for skewed data
  if (skew) {
    
    # score distances using adjusted outlyingness
    sd <- outl2
    
    cutoff.sd <- coOD(sd, skew=TRUE)$co.od
  }
  flag.sd <- (sd<=cutoff.sd) # FALSE if outlying in PCA subspace
  flag.od <- (od<=cutoff.od) # FALSE if outlying for OD
  flag.all <- (flag.od*flag.sd)==1 # FALSE if outlier (all 3 types)
  
  return(list(loadings=P, eigenvalues=l, scores=Tn, center=mu, D=D, k=k, H0=H0,
              H1=H1, P1=P1, index=index, H2=H2, P2=P2, H3=H3, alpha=alpha, h=h,
              sd=sd, od=od, cutoff.sd=cutoff.sd, cutoff.od=cutoff.od,
              flag.sd=flag.sd, flag.od=flag.od, flag.all=flag.all))
}