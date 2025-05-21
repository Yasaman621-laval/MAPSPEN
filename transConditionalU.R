`transConditionalU` <- function(summaryZ, Nvec, JointBmatrix, Zcov, SDvec= NULL, Zscale = 1, sigma2K_allAlpha_List, es_alphaMatrix, plinkLD, weight = 1.0){
 

  Nq = length(Nvec)


  if(is.null(SDvec)){
    summaryBetas = matrix(0, nrow(summaryZ), Nq)
    rownames(summaryBetas) =  rownames(summaryZ)
    SDvec = matrix(0, nrow(summaryZ), Nq)
    for(ii in 1:Nq){
        summaryBetas[,ii] = summaryZ[,ii]/sqrt(Nvec[ii])
        SDvec[,ii] = 1/sqrt(Nvec[ii])
    }
  }else{
    summaryBetas = summaryZ
  }
  
  Betaindex = c(1:nrow(summaryBetas))-1
  SNPnames = rownames(summaryBetas)
  
  ldJ = PlinkLD_transform_addR(plinkLD, SNPnames)
  rm(plinkLD)

  JidMatrix = matrix(,nrow(ldJ),2)
  mat1 = match(ldJ[,1],SNPnames)
  mat2 = match(ldJ[,2],SNPnames)

  JidMatrix[,1] = Betaindex[mat1]
  JidMatrix[,2] = Betaindex[mat2]

  ldJ[,1] = JidMatrix[,1]
  ldJ[,2] = JidMatrix[,2]
  
  ldJIndexMatrix = ldJ[,1:2]
  ldJRMatrix = ldJ[,3:4]

  od = order(JidMatrix[,1],JidMatrix[,2], decreasing=F)
  ldJ = ldJ[od,]
  
  wind = which(! Betaindex %in% ldJ[,1])

  IndJ = -1

  if(length(wind) > 0){
    IndJ = Betaindex[wind]
  }

  Counts = table(ldJ[,1])
  NumSNP = length(Counts)
  
  IndexS = c(0,cumsum(Counts)[-NumSNP])
  IndexE = cumsum(Counts)-1

  IndexMatrix = matrix(,NumSNP,3)
  IndexMatrix[,1] = as.numeric(names(Counts))

  nrow_IndexMatrix = NumSNP
  ncol_IndexMatrix = ncol(IndexMatrix)


  IndexMatrix[,2] = IndexS
  IndexMatrix[,3] = IndexE


  ldvec = ldJ[,3:4]
  ldJ = ldJ[,2]

  length_ldJ = length(ldJ)



  P = nrow(summaryBetas)
  Q = ncol(summaryBetas)
  
  nrow_JointBmatrix = nrow(JointBmatrix)
  ncol_JointBmatrix = ncol(JointBmatrix)
  
  Zcov_Inv = solve(Zcov)
  nrow_Zcov_Inv = nrow(Zcov_Inv)
  ncol_Zcov_Inv = ncol(Zcov_Inv)
  
  det_Zcov = det(Zcov)
  sqrt_det_Zcov = sqrt(det_Zcov)
  
  
  numAlpha = 2^Q
  num_element = Q*Q
  SharedPattern = permutations(n=2,r=Q,v=c(0,1),repeats.allowed=T)
  numI = apply(SharedPattern,1,sum)
  Nmatrix = matrix(0,Q,Q)
  diag(Nmatrix) = Nvec
  A0 = sqrt(Nmatrix) %*% Zcov_Inv %*% sqrt(Nmatrix)
  A0_inv = solve(A0)
  AllAn_i = matrix(0, numAlpha,num_element)


  sqrt_det_Ai = rep(1,numAlpha)
  sqrt_inv_det_An_i = rep(1,numAlpha)


  for(aIndex in 1:numAlpha){
    if(numI[aIndex]!=0){
      wI = which(SharedPattern[aIndex,]==1)
      Sigma_i = sigma2K_allAlpha_List[[aIndex]]
      Sigma_i = Sigma_i[wI,wI]
      A_i = solve(Sigma_i)
      A0_i = A0[wI,wI]
      An_i = (A0_i + A_i)
      An_i_inv = solve(An_i)

      k = 1
      for(a1 in 1:numI[aIndex]){
        for(a2 in 1:numI[aIndex]){
          AllAn_i[aIndex,k] = An_i_inv[a1,a2]
          k = k + 1
        }
      }

      sqrt_det_Ai[aIndex] = sqrt(det(A_i))
      sqrt_inv_det_An_i[aIndex] = 1/sqrt(det(An_i))
    }

  }

  nrow_es_alphaMatrix = nrow(es_alphaMatrix)
  ncol_es_alphaMatrix = ncol(es_alphaMatrix)


  nrow_SharedPattern = nrow(SharedPattern)
  ncol_SharedPattern = ncol(SharedPattern)

  DensityU = matrix(0, nrow_JointBmatrix, numAlpha*P)
  
  b0Matrix =  matrix(0, nrow(summaryBetas), Q)
  esjointbetas = matrix(0, nrow(summaryBetas), Q)



  dims = rep(0,21)
  dims[1] = nrow(ldvec)
  dims[2] = P
  dims[3] = nrow_JointBmatrix
  dims[4] = ncol_JointBmatrix
  dims[5] = nrow_IndexMatrix
  dims[6] = ncol_IndexMatrix
  dims[7] = length(wind)
  dims[8] = Zscale
  dims[11] = Q
  dims[12] = numAlpha
  dims[13] = num_element
  dims[14] = nrow_es_alphaMatrix
  dims[15] = ncol_es_alphaMatrix
  dims[18] = nrow_SharedPattern
  dims[19] = ncol_SharedPattern
  dims[20] = nrow(DensityU)
  dims[21] = ncol(DensityU)


  numvec = rep(0,1)
  numvec[1] = 1/sqrt_det_Zcov

  postEb1 = matrix(0,nrow_JointBmatrix, numAlpha*P)
  postEb2 = matrix(0,nrow_JointBmatrix, numAlpha*P)


  Z = .C("transConditionalU", as.double(t(summaryBetas)),as.integer(ldJ), as.integer(dims), as.integer(t(IndexMatrix)), as.integer(IndJ),
      as.double(t(ldvec)), as.double(t(SDvec)),DensityU = as.double(t(DensityU)),
      as.double(t(JointBmatrix)),as.double(t(AllAn_i)), as.double(Nvec), as.double(sqrt_det_Zcov),as.double(sqrt_det_Ai),as.double(sqrt_inv_det_An_i),as.double(t(A0)), as.double(t(A0_inv)), as.double(t(es_alphaMatrix)),as.integer(t(SharedPattern)),as.integer(numI), as.double(weight), as.double(Zcov_Inv), b0Matrix = as.double(t(b0Matrix)), esjointbetas = as.double(t(esjointbetas)), postEb1 = as.double(t(postEb1)), postEb2 = as.double(t(postEb2)), as.double(numvec),PACKAGE="SummaryLasso")
    
    

  DensityU = matrix(Z$DensityU, nrow = nrow_JointBmatrix, ncol = numAlpha*P, byrow = TRUE)
  b0Matrix = matrix(Z$b0Matrix, nrow = nrow(summaryBetas), ncol = Q, byrow = TRUE)
  esjointbetas = matrix(Z$esjointbetas, nrow = nrow(summaryBetas), ncol = Q, byrow = TRUE)
  postEb1 = matrix(Z$postEb1, nrow = nrow_JointBmatrix, ncol = numAlpha*P, byrow = TRUE)
  postEb2 = matrix(Z$postEb2, nrow = nrow_JointBmatrix, ncol = numAlpha*P, byrow = TRUE)


  
  ll = list(DensityU = DensityU, b0Matrix = b0Matrix, esjointbetas = esjointbetas, postEb1 = postEb1, postEb2 = postEb2)
  return(ll)
  
  
   
}

