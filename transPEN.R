`transPEN` <- function(summaryZ, summaryBetas = NULL, SSvec = NULL, Nvec, plinkLDList, NumIter = 100, breaking = 1,
  numChrs=22, ChrIndexBeta = 0, Init_summaryBetas= NULL, Zscale = 1, RupperVal = NULL,
  tuningMatrix = NULL, penalty=c("Lasso","LassowII","LassowIII","transw","mixLOG"),
  taufactor = c(1/25, 1, 10), llim_length =2, subtuning = 50, Lambda_limit = c(0.5,0.9), 
  Lenlam_singleTrait = 200, dfMax = NULL, IniBeta = 0, inverseTuning = 0, outputAll = 0, warmStart = 1 , orderseq = 1, basew = NULL, tau = 0.0, alpha = 1, singleStart = 1, weightN = 0){
 

  #if(Zscale!=1){stop("Tuning values set-up for multiple traits analysis requires Zscale=1.")}

  if(inverseTuning==1){
    orderD = "FALSE"
  }else{
    orderD = "TRUE"
  }
  Nq = length(Nvec)
  
  
  
  if(is.null(summaryBetas)){
    summaryBetas = matrix(0, nrow(summaryZ), Nq)
    SDvec = matrix(0, nrow(summaryZ), Nq)
  
    for(ii in 1:Nq){
      summaryBetas[,ii] = summaryZ[,ii]/sqrt(Nvec[ii])
      SDvec[,ii] = 1/sqrt(Nvec[ii])
    }
    rownames(summaryBetas) = rownames(summaryZ)
    SSvec = matrix(0,1,1)
  }else{
    SDvec = matrix(0,1,1)
  }
  
  if(!is.null(Init_summaryBetas)){
    IniBeta = 1
  }else{
    Init_summaryBetas = matrix(0, nrow(summaryBetas), Nq)
  }
   

  
  sumdiffbetas = apply(summaryBetas,1,abs_diff)

  

  if(is.null(dfMax)){
    dfMax = ceiling(0.7*nrow(summaryBetas))
  }
  
  autoTuning = 0
  if(is.null(tuningMatrix)){
    autoTuning = 1
    tauvec = tau
    Lambda_limit = quantile(abs(summaryZ[,1]), Lambda_limit)

    tuningMatrix = Tuning_setup_group_only(tauvec, subtuning, Lambda_limit, llim_length, Nq)
  }            
     

  inv_summaryBetas = 0

  
  Betaindex = c(1:nrow(summaryBetas))-1
  SNPnames = rownames(summaryBetas)


  q_nrowIndexMatrix = q_ldJ = q_IndJ = q_ldvec = c()
  
  for(iiIndex in 1:Nq){
    ldJ = PlinkLD_transform(plinkLDList[[iiIndex]], SNPnames)

    JidMatrix = matrix(,nrow(ldJ),2)
    mat1 = match(ldJ[,1],SNPnames)
    mat2 = match(ldJ[,2],SNPnames)

    JidMatrix[,1] = Betaindex[mat1]
    JidMatrix[,2] = Betaindex[mat2]

    ldJ[,1] = JidMatrix[,1]
    ldJ[,2] = JidMatrix[,2]

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

    ldvec = ldJ[,3]
    ldJ = ldJ[,2]

    if(iiIndex==1){
      All_IndexMatrix = IndexMatrix
      All_ldJ = ldJ
      All_IndJ = IndJ
      All_ldvec = ldvec
      All_length_wind = length(wind)
    }else{
      All_IndexMatrix = rbind(All_IndexMatrix, IndexMatrix)
      All_ldJ = c(All_ldJ,ldJ)
      All_IndJ = c(All_IndJ,IndJ)
      All_ldvec = c(All_ldvec,ldvec)
      All_length_wind = c(All_length_wind,length(wind))
    }
    

    q_nrowIndexMatrix = c(q_nrowIndexMatrix,NumSNP)
    q_ldJ = c(q_ldJ, length(ldJ))
    q_IndJ = c(q_IndJ, length(IndJ))
    q_ldvec = c(q_ldvec, length(ldvec))
    rm(ldJ,IndJ,ldvec,wind)
  }

  if(is.null(RupperVal)){
    RupperVal = ceiling(max(abs(summaryBetas),na.rm=T)*50)
  }

  P = nrow(summaryBetas)
  Q = ncol(summaryBetas)
  
  if(nrow(tuningMatrix)>1){
    NumTuning = nrow(tuningMatrix)
  }
  
  if(warmStart==1){
    if(Nq>1){
      if(orderseq==1){
        if(penalty=="Lasso"){
          id = order(tuningMatrix[,1],1/tuningMatrix[,2],decreasing=orderD)
        }else{
          id = order(tuningMatrix[,1],tuningMatrix[,2],decreasing=orderD)
        }
        useI = 1
      }
      if(orderseq==2){
        if(penalty=="Lasso"){
          id = order(1/tuningMatrix[,2],tuningMatrix[,1],decreasing=orderD)
        }else{
          id = order(tuningMatrix[,2],tuningMatrix[,1],decreasing=orderD)
        }
        useI = 2
      }
      tuningMatrix = tuningMatrix[id,]
      groupvec = apply(tuningMatrix[,c(useI),drop=F],1,paste0,collapse="")
      Tgroupvec = table(groupvec)
      Ugroup = unique(groupvec)
      mat = match(Ugroup, names(Tgroupvec))
      Tgroupvec = Tgroupvec[mat]
      if(length(Ugroup)==1){
        StartVec = c(1)
      }else{
        StartVec = c(1, (cumsum(Tgroupvec)[1:(length(cumsum(Tgroupvec))-1)]+1))
      }
    }else{
      if(orderseq==1){
        id = order(tuningMatrix[,1],decreasing=orderD)
      }
      if(orderseq==2){
        id = order(tuningMatrix[,3],decreasing=orderD)
      }
      tuningMatrix = tuningMatrix[id,]
      StartVec = 1
    }
    tuningMatrix = cbind(tuningMatrix,0)
    tuningMatrix[StartVec,ncol(tuningMatrix)] = 1
  }else{
    tuningMatrix = cbind(tuningMatrix,0)
  }
  if(singleStart==1){
    tuningMatrix[,ncol(tuningMatrix)] = 0
    tuningMatrix[1,ncol(tuningMatrix)] = 1
  }

  nrow_All_IndexMatrix = nrow(All_IndexMatrix)
  ncol_All_IndexMatrix = ncol(All_IndexMatrix)

  ncolBetaMatrix = P*Q
  dims = rep(0,20)
  dims[1] = NumTuning
  dims[2] = P
  dims[3] = NumIter
  dims[4] = breaking
  dims[5] = nrow_All_IndexMatrix
  dims[6] = ncol_All_IndexMatrix
  dims[7] = 0  #NumInd
  dims[8] = Zscale
  dims[9] = nrow(tuningMatrix)
  dims[10] = ncol(tuningMatrix)
  dims[11] = Q
  dims[12] = ncolBetaMatrix
  dims[13] = outputAll
  dims[14] = dfMax
  dims[15] = warmStart
  dims[16] = IniBeta
  dims[17] = nrow(SDvec)
  dims[18] = ncol(SDvec)
  dims[19] = nrow(SSvec)
  dims[20] = ncol(SSvec)
  
  if(is.null(basew)){
    basew = max(sumdiffbetas,na.rm=T)
  }
  doublevec = c(RupperVal, basew, tau, alpha, weightN)

  Numitervec = rep(0, NumTuning)
  Dfq1 = rep(0, NumTuning)
  BetaMatrix = matrix(0, NumTuning, ncolBetaMatrix)
  
  
  

  Z = .C("transPEN", as.double(t(summaryBetas)),as.integer(All_ldJ), as.integer(dims), Numitervec=as.integer(Numitervec),
   as.integer(t(All_IndexMatrix)), as.integer(All_IndJ),
   as.double(All_ldvec), as.double(inv_summaryBetas), as.integer(ChrIndexBeta), as.double(doublevec),
   as.double(t(Init_summaryBetas)), as.double(t(SDvec)),
   as.double(t(tuningMatrix)), BetaMatrix = as.double(t(BetaMatrix)),
   penalty, as.integer(q_nrowIndexMatrix), as.integer(q_ldJ),as.integer(q_IndJ),as.integer(q_ldvec), as.integer(All_length_wind), Dfq1 = as.integer(Dfq1), as.double(t(SSvec)),PACKAGE="SummaryLasso")




  BetaMatrix = matrix(Z$BetaMatrix, nrow = NumTuning, ncol = ncolBetaMatrix, byrow = TRUE)
  #tuningMatrix = tuningMatrix[,-ncol(tuningMatrix)]
  
  #if(autoTuning){
   # tuningMatrix[c(1:Lenlam_singleTrait), c(2:3)] = NA
  #} 
  colnames(BetaMatrix) = paste0(rep(SNPnames, times = Q),".trait", rep(c(1:Q), each = P))
   
  Numitervec = Z$Numitervec
  Dfq1 = Z$Dfq1
   
  if(outputAll==0){
    convergeIndex = which(Numitervec > 0)
    Numitervec = Numitervec[convergeIndex]
    BetaMatrix = BetaMatrix[convergeIndex,]
    tuningMatrix = tuningMatrix[convergeIndex,]
  }
   
  ll = list(BetaMatrix = BetaMatrix, Numitervec = Numitervec, tuningMatrix = tuningMatrix,Dfq1=Dfq1)
  return(ll)
   
}

