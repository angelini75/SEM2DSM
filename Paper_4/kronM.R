# RHO=matrix(c(1,10,100,
#            10, 1,10,
#            100,10,1), 3,3)
# dimnames(RHO) <- list(letters[10:12], letters[10:12]) 
# 
# RHO
# #     j  k   l
# # j   1 10 100
# # k  10  1  10
# # l 100 10   1
# 
# RHO.I=matrix(c(1,0,0,
#             0,1,0,
#             0,0,1), 3,3)
# dimnames(RHO.I) <- list(letters[10:12], letters[10:12]) 
# 
# RHO.I
# #   j k l
# # j 1 0 0
# # k 0 1 0
# # l 0 0 1
# 
# SIGMA0= matrix(rep(7:8),4,4)
# diag(SIGMA0) <- 1:4
# dimnames(SIGMA0) <- list(LETTERS[1:4], LETTERS[1:4]) 
# 
# SIGMA0
# #   A B C D
# # A 1 7 7 7
# # B 8 2 8 8
# # C 7 7 3 7
# # D 8 8 8 4
# 

kronM <- function (RHO, RHO.I, SIGMA0, sp = 1:9)#, make.dimnames = TRUE) 
{
  RHO <- as.array(RHO)
  RHO.I <- as.array(RHO.I)
  SIGMA0 <- as.array(SIGMA0)
  dnx1 <- dimnames(RHO)
  dnx2 <- dimnames(RHO.I)
  dny <- dimnames(SIGMA0)
  dRHO <- dim(RHO)
  dRHO.I <- dim(RHO.I)
  dSIGMA0 <- dim(SIGMA0)
  #ld <- length(dX) - length(dSIGMA0)
  opobjRHO <- outer(RHO, SIGMA0)
  opobjRHO.I <- outer(RHO.I, SIGMA0)
  opobj <- opobjRHO.I
  opobj[,,sp,sp] <- opobjRHO[,,sp,sp]
  
  dp <- as.vector(t(matrix(1L:(2 * length(dRHO)), ncol = 2)[, 
                                                          2:1]))
  opobj <- aperm(a = opobj, dp)
  dim(opobj) <- dRHO * dSIGMA0
  # if (make.dimnames && !(is.null(dnx) && is.null(dny))) {
  #   if (is.null(dnx)) 
  #     dnx <- vector("list", length(dX))
  #   else if (ld < 0L) 
  #     dnx <- c(dnx, vector("list", -ld))
  #   tmp <- which(sapply(dnx, is.null))
  #   dnx[tmp] <- lapply(tmp, function(i) rep.int("", dX[i]))
  #   if (is.null(dny)) 
  #     dny <- vector("list", length(dSIGMA0))
  #   else if (ld > 0) 
  #     dny <- c(dny, vector("list", ld))
  #   tmp <- which(sapply(dny, is.null))
  #   dny[tmp] <- lapply(tmp, function(i) rep.int("", dSIGMA0[i]))
  #   k <- length(dim(opobj))
  #   dno <- vector("list", k)
  #   for (i in 1L:k) {
  #     tmp <- outer(dnx[[i]], dny[[i]], FUN = "paste", sep = ":")
  #     dno[[i]] <- as.vector(t(tmp))
  #   }
  #   dimnames(opobj) <- dno
  # }
  opobj
}