PlinkLD_transform_addR <- function(plinkLD, Allkeepsnps){

  greR = grep("R.",colnames(plinkLD),fixed=T)
  namesR = colnames(plinkLD)[greR]
  namesR = sort(namesR,decreasing = FALSE)

  wkeep = which( plinkLD$SNP_A  %in% Allkeepsnps & plinkLD$SNP_B %in% Allkeepsnps)
  plinkLD = plinkLD[wkeep,]

  if(length(which(is.na(plinkLD)))>0){stop("")}
	
  ldJ = plinkLD[,c("SNP_B","SNP_A",namesR)]
  names(ldJ) = c("SNP_A","SNP_B",namesR)
	  
  ldJ = rbind(plinkLD[,c("SNP_A","SNP_B",namesR)],ldJ)

  return(ldJ)
}
