# SNP testing pupation and eclosion
# independent and normalised by HPD

##################
# Data cleaning and set up
##################
library(gtools)
library(limma)

load("../Data/alldiets_mean_percentage.RData")

alldiets_mean_pupation_percentage_HPD = 
  alldiets_mean_pupation_percentage - alldiets_mean_pupation_percentage[,"HPD"]
alldiets_mean_pupation_percentage_HPD = alldiets_mean_pupation_percentage_HPD[,c("HSD","HFD","HFDlard","WD","HSt")]
# taking the difference of the percentage compared to protein
# so a positive value means it's higher in the diet compared to protein

alldiets_mean_eclosion_percentage_HPD = 
  alldiets_mean_eclosion_percentage - alldiets_mean_eclosion_percentage[,"HPD"]
alldiets_mean_eclosion_percentage_HPD = alldiets_mean_eclosion_percentage_HPD[,c("HSD","HFD","HFDlard","WD","HSt")]

save(alldiets_mean_pupation_percentage_HPD,alldiets_mean_eclosion_percentage_HPD,
     file = "alldiets_mean_percentage_HPD.RData")


for (type in c("", "_HPD")) {
  
  if (type == "") {
    
    load("../Data/alldiets_mean_percentage.RData")
    
  }
  
  if (type == "_HPD") {
    
    load("alldiets_mean_percentage_HPD.RData")
    alldiets_mean_pupation_percentage = alldiets_mean_pupation_percentage_HPD
    alldiets_mean_eclosion_percentage = alldiets_mean_eclosion_percentage_HPD
    
  }
  
  
  # pupation and eclosion matrices have the same dimension
  alldiets_mean_pupation_percentage
  alldiets_mean_eclosion_percentage
  lines = mixedsort(rownames(alldiets_mean_pupation_percentage))
  
  dgrp_raw = as.matrix(read.table("../DGRP_Data/dgrp2.tgeno",header=TRUE))
  rownames(dgrp_raw) <- dgrp_raw[,"id"]
  dgrp_raw <- dgrp_raw[mixedsort(rownames(dgrp_raw)),]
  save(dgrp_raw,file="../DGRP_Data/dgrp_raw.RData")
  snps = rownames(dgrp_raw)
  
  dgrp_annot = as.matrix(read.table("../DGRP_Data/dgrp.fb557.annot.txt",header=FALSE))
  colnames(dgrp_annot) <- c("id","SNP","siteclass","tf_binding")
  rownames(dgrp_annot) <- dgrp_annot[,"id"]
  dgrp_annot <- dgrp_annot[snps,]
  save(dgrp_annot,file="../DGRP_Data/dgrp_annot.RData")
  
  dgrp_geno = dgrp_raw[,lines]
  save(dgrp_geno,file="../DGRP_Data/dgrp_geno.RData")
  
  wolbachia_raw = read.table("../DGRP_Data/wolbachia.tsv",header = FALSE, row.names = 1)
  wolbachia = wolbachia_raw[,1]
  names(wolbachia) = rownames(wolbachia_raw)
  wolb <- wolbachia[lines]
  save(wolb,file="../DGRP_Data/wolbachia.RData")
  
  # line relationship matrix
  rel_raw = read.delim("../DGRP_Data/freeze2.common.rel.mat.txt", row.names = 1, 
                       header = TRUE)
  # extract the first 10 PCs
  rel_pc = princomp(rel_raw)
  rel = rel_pc$scores[,1:10]
  rownames(rel) <- gsub(" .*", "", rownames(rel))
  table(rownames(rel) %in% lines)
  table(lines %in% rownames(rel))
  rel_lines = rel[lines,]
  saveRDS(rel_lines, file = "../DGRP_Data/rel_lines.Rds")
  
  inversion_raw = read.table("../DGRP_Data/inversion.txt", header = TRUE, 
                             row.names = 1, stringsAsFactors = TRUE)
  
  dim(inversion_raw)
  summary(inversion_raw)
  rownames(inversion_raw) <- gsub("DGRP", "line", rownames(inversion_raw))
  
  # get the line names
  inversion <- inversion_raw[lines,]
  inversion <- apply(inversion, 2, function(x){
    ifelse(x == "ST","ST","Other")
  })
  inversion <- as.data.frame(inversion)
  
  minLinesInversion = sort(apply(inversion, 2, function(x) min(table(factor(x, levels = c("ST","Other"))))))
  
  minLinesInversion
  # there are only 5 inversion statuses with 5 or more lines exhibiting
  # these go ahead as covariates in the MANOVA test
  
  inversionSelectedLines = names(which(minLinesInversion >= 5))
  inversionSelectedLines
  
  inversion_selected = inversion[,inversionSelectedLines]
  
  dim(inversion_selected)
  saveRDS(inversion_selected, file = "../DGRP_Data/dgrp_inversions.Rds")
  
  load("../Data/processedPhenotypeData.RData")
  batch = setNames(HPD_processed$clean[lines, "Batch"], lines)
  
  identical(rownames(alldiets_mean_pupation_percentage),colnames(dgrp_geno))
  identical(names(wolb),colnames(dgrp_geno))
  identical(names(wolb),rownames(inversions))
  
  minAlleles = 5
  
  inversions_sub = inversions[,apply(inversions, 2, function(x)min(table(x)) >= minAlleles)]
  wilcox_data = cbind(wolb, inversions_sub, rel, factor(batch))
  
  alldiets_mean_pupation_percentage_res = apply(alldiets_mean_pupation_percentage, 2, function(y) {
    lm(y ~ ., data = wilcox_data)$residuals[lines]
  })
  
  alldiets_mean_eclosion_percentage_res = apply(alldiets_mean_eclosion_percentage, 2, function(y) {
    lm(y ~ ., data = wilcox_data)$residuals[lines]
  })
  
  ##################
  # manova testing
  ##################
  
  # check these two are TRUE
  identical(rownames(alldiets_mean_pupation_percentage),colnames(dgrp_geno))
  identical(names(wolb),colnames(dgrp_geno))
  
  manova_pvalPupation_batch = rep(NA,nrow(dgrp_geno))
  names(manova_pvalPupation_batch) <- rownames(dgrp_geno)
  manova_pvalEclosion_batch = rep(NA,nrow(dgrp_geno))
  names(manova_pvalEclosion_batch) <- rownames(dgrp_geno)
  
  manova_pvalPupation_wolb = rep(NA,nrow(dgrp_geno))
  names(manova_pvalPupation_wolb) <- rownames(dgrp_geno)
  manova_pvalEclosion_wolb = rep(NA,nrow(dgrp_geno))
  names(manova_pvalEclosion_wolb) <- rownames(dgrp_geno)
  
  manova_coefPupation_batch = matrix(NA, nrow = nrow(dgrp_geno),
                                     ncol = ncol(alldiets_mean_pupation_percentage),
                                     dimnames = list(
                                       rownames(dgrp_geno),
                                       colnames(alldiets_mean_pupation_percentage)
                                     ))
  
  manova_coefPupation_wolb <- manova_coefPupation_batch
  manova_coefEclosion_batch <- manova_coefPupation_batch
  manova_coefEclosion_wolb <- manova_coefPupation_batch
  
  for (i in 1:nrow(dgrp_geno)) {
    if (i%%1000==0) print(i)
    x = dgrp_geno[i,]
    if (min(sum(x=="0"),sum(x=="2"))<minAlleles) {
      next
    } else {
      inds = x=="0"|x=="2"
      xs = x[inds]
      ysPupation = alldiets_mean_pupation_percentage[inds,]
      ysEclosion = alldiets_mean_eclosion_percentage[inds,]
      wolbs = wolb[inds]
      inversions_inds = inversions[inds,]
      rels = rel[inds,]
      batchs = factor(batch[inds])
      inversions_inds <- inversions_inds[,apply(inversions_inds,
                                                2,
                                                function(x)min(table(x)) >= minAlleles)]
      all_data_batch = cbind(xs, wolbs, inversions_inds, rels, batchs)
      all_data_wolb = data.frame(cbind(xs, wolbs))
      fitPupation_batch = 
        tryCatch(manova(ysPupation~., data = all_data_batch), 
                 error=function(e) return(NA))
      fitEclosion_batch = 
        tryCatch(manova(ysEclosion~., data = all_data_batch), 
                 error=function(e) return(NA))
      fitPupation_wolb = 
        tryCatch(manova(ysPupation~., data = all_data_wolb), 
                 error=function(e) return(NA))
      fitEclosion_wolb = 
        tryCatch(manova(ysEclosion~., data = all_data_wolb), 
                 error=function(e) return(NA))
      if (!is.na(fitPupation_batch)) {
        p = tryCatch(summary(fitPupation_batch, test="Wilks", tol = 0)$stats[1,"Pr(>F)"],
                     error=function(e) return(NA))
        if (is.null(p)|anyNA(p[1])) {
          p <- NA
        } else {
          manova_pvalPupation_batch[i] = p
          
          manova_coefPupation_batch[i,] <- fitPupation_batch$coefficients["xs2",]
          
        }
      }
      
      if (!is.na(fitEclosion_batch)) {
        
        p = tryCatch(summary(fitEclosion_batch, test="Wilks", tol = 0)$stats[1,"Pr(>F)"],
                     error=function(e) return(NA))
        
        if (is.null(p)|anyNA(p[1])) {
          p <- NA
        } else {
          manova_pvalEclosion_batch[i] = p
          manova_coefEclosion_batch[i,] <- fitEclosion_batch$coefficients["xs2",]
        }
      }
      
      if (!is.na(fitPupation_wolb)) {
        p = tryCatch(summary(fitPupation_wolb, test="Wilks", tol = 0)$stats[1,"Pr(>F)"],
                     error=function(e) return(NA))
        
        if (is.null(p)|anyNA(p[1])) {
          p <- NA
        } else {
          manova_pvalPupation_wolb[i] = p
          
          manova_coefPupation_wolb[i,] <- fitPupation_wolb$coefficients["xs2",]
        }
      }
      
      if (!is.na(fitEclosion_wolb)) {
        p = tryCatch(summary(fitEclosion_wolb, test="Wilks", tol = 0)$stats[1,"Pr(>F)"],
                     error=function(e) return(NA))
        
        if (is.null(p)|anyNA(p[1])) {
          p <- NA
        } else {
          manova_pvalEclosion_wolb[i] = p
          manova_coefEclosion_wolb[i,] <- fitEclosion_wolb$coefficients["xs2",]
        }
      }
    }
  }
  
  numReference = apply(dgrp_geno == "0", 1, sum, na.rm = TRUE)
  numAlternate = apply(dgrp_geno == "2", 1, sum, na.rm = TRUE)
  numList = list(numReference = numReference, numAlternate = numAlternate)
  saveRDS(numList, file = "../numList.Rds")
  
  out = list(
    manova_pvalPupation_batch = manova_pvalPupation_batch,
    manova_coefPupation_batch = manova_coefPupation_batch,
    
    manova_pvalPupation_wolb = manova_pvalPupation_wolb,
    manova_coefPupation_wolb = manova_coefPupation_wolb,
    
    manova_pvalEclosion_batch = manova_pvalEclosion_batch,
    manova_coefEclosion_batch = manova_coefEclosion_batch,
    
    manova_pvalEclosion_wolb = manova_pvalEclosion_wolb,
    manova_coefEclosion_wolb = manova_coefEclosion_wolb,
    
    numReference = numReference,
    numAlternate = numAlternate
    
  )
  saveRDS(out, file = paste0("../manova_pvals_coefs", type, ".Rds"))
  
}