library(reshape)
library(ggplot2)
library(viridis)
library(gtools)
library(ggthemes)
library(dplyr)
library(xlsx)

lineConversion = read.xlsx("../Data/DGRP line information.xlsx",
                           sheetIndex = 1,header=FALSE)
colnames(lineConversion) <- c("Essi_line","Stock_number","line")
# "line" is the line that DGRP uses and I use
lineConversion[,"line"] <- gsub("DGRP-","line_",lineConversion[,"line"])

diets = c("HPD","HSD","HFD","HFDlard","WD","HSt")
phenotypeFile = "../Data/DGRP data 200917_FOR SHILA_ver4_SG.xlsx"
# note the change of sheet name to HFD lard to HFDlard

processRawPhenotypeData = function(raw) {
  # fix column ids and add the "line" column to the dataset
  clean = raw
  colnames(clean)[1:3] <- c("Batch","Essi_line","Stock_number")
  rownames(clean) <- lineConversion[match(clean[,"Stock_number"],
                                          lineConversion[,"Stock_number"]),"line"]
  
  colnames(clean)[4:31] <- apply(expand.grid(paste0("Rep",c("A","B","C","D")),
                                             paste0("Day",as.character(5:11))),
                                 1,paste0,collapse="_")
  
  colnames(clean)[32:39] <- apply(expand.grid(paste0("Rep",c("A","B","C","D")),
                                              c("NonEclosed","Eclosed")),
                                  1,paste0,collapse="_")
  
  # if data is missing for day 11, it is assumed to be zero
  # ensure that 'missing' reps remain missing though
  
  clean[,grep("Day11",colnames(clean))][is.na(clean[,grep("Day11",colnames(clean))])] <- 0
  #image(is.na(clean))
  
  # put day 11 cells that should be missing back in
  
  rowcolnames = t(apply(which(is.na(clean[,grep("Day5",colnames(clean))]),arr=TRUE),1,
        function(x)c(rownames(clean)[x[1]],paste0("Rep",c("A","B","C","D"),"_Day11")[x[2]])))
  for (i in 1:nrow(rowcolnames)) {
    clean[rowcolnames[i,1],rowcolnames[i,2]] <- NA
  }
  
  # calculate number (out of the 30 eggs) that pupated
  # this should be a lines x 4 matrix (for A B C D)
  pupation_count = sapply(c("A","B","C","D"),
                          function(l) rowSums(clean[,grepl(paste0("Rep",l,"_Day"),
                                                           colnames(clean))]))
  pupation_percentage = pupation_count/30
  mean_pupation_percentage = rowMeans(pupation_percentage,na.rm=TRUE)
  
  # clean up the eclosion set
  for (j in rownames(clean)) {
    print(j)
    for (k in colnames(clean)[grepl("Eclosed",colnames(clean))]) {
      if (!is.na(clean[j,k])) next
      repl = substring(k,4,4)
      type = gsub(".+_","",k)
      pupated = pupation_count[j,repl]
      othertype = ifelse(type=="Eclosed","NonEclosed","Eclosed")
      othercount = clean[j,paste0("Rep",repl,"_",othertype)]
      clean[j,k] <- pupated-othercount
    }
  }
  
  eclosed_count = as.matrix(clean[,grepl("_Eclosed",colnames(clean))])
  colnames(eclosed_count) <- substring(colnames(eclosed_count),4,4)
  eclosed_percentage = eclosed_count/pupation_count
  mean_eclosed_percentage = rowMeans(eclosed_percentage,na.rm=TRUE)
  
  m = melt(pupation_percentage)
  colnames(m) <- c("line","Rep","Pupation_percentage")
  fit_rep_pup = lm(Pupation_percentage ~ line + Rep,data=m)
  anova(fit_rep_pup)
  # observe significant line effect but not batch effect
  
  boxplot(mean_pupation_percentage~clean[,"Batch"],
          main = "Mean pupation percentage")
  fit_batch_pup = lm(mean_pupation_percentage~factor(clean[,"Batch"]))
  anova(fit_batch_pup)
  # no single batch is significantly (P<0.05) higher than the rest
  
  me = melt(eclosed_percentage)
  colnames(me) <- c("line","Rep","Eclosed_percentage")
  fit_rep_ecl = lm(Eclosed_percentage ~ line + Rep,data=me)
  anova(fit_rep_ecl)
  # observe significant line effect but not batch effect
  
  boxplot(mean_eclosed_percentage~clean[,"Batch"],
          main = "Mean eclosion percentage")
  fit_batch_ecl = lm(mean_eclosed_percentage~factor(clean[,"Batch"]))
  anova(fit_batch_ecl)
  # no single batch is significantly (P<0.05) higher than the rest
  
  return(list(clean=clean,
              pupation_count=pupation_count,
              pupation_percentage = pupation_percentage,
              mean_pupation_percentage = mean_pupation_percentage,
              fit_rep_pup = fit_rep_pup,
              fit_batch_pup = fit_batch_pup,
              fit_rep_ecl = fit_rep_ecl,
              fit_batch_ecl = fit_batch_ecl,
              eclosed_count = eclosed_count,
              eclosed_percentage = eclosed_percentage,
              mean_eclosed_percentage = mean_eclosed_percentage
              ))
}

for (i in 1:length(diets)) {
  print(i)
  assign(paste0(diets[i],"_raw"),read.xlsx(phenotypeFile,sheetName = diets[i],
                                         startRow=4,endRow=200,colIndex = 1:39))
}
save(list=paste0(diets,"_raw"),file="../Data/rawPhenotypeData.RData")

for (i in 1:length(diets)) {
  print(i)
  assign(paste0(diets[i],"_processed"),processRawPhenotypeData(get(paste0(diets[i],"_raw"))))
}
save(list=paste0(diets,"_processed"),file="../Data/processedPhenotypeData.RData")

# relationship between pupation and eclosion percentages by batch
pdf(file="../Figures/boxplots_pupation_eclosion_batch.pdf",height=6,width=12)
par(mfcol=c(2,6))
par(mar=c(5,4,2,1))
for (i in 1:length(diets)) {
  print(i)
  processed = get(paste0(diets[i],"_processed"))
  boxplot(processed$mean_pupation_percentage~processed$clean[,"Batch"],
          main = diets[i],ylab="Mean pupation percentage",col=i+1,xlab="Batch")
  boxplot(processed$mean_eclosed_percentage~processed$clean[,"Batch"],
          ylab = "Mean eclosion percentage",col=i+1,xlab="Batch")
}
dev.off()

alldiets_mean_pupation_percentage = sapply(
  diets,function(diet)get(paste0(diet,"_processed"))$mean_pupation_percentage)
alldiets_mean_eclosion_percentage = sapply(
  diets,function(diet)get(paste0(diet,"_processed"))$mean_eclosed_percentage)

par(mfrow=c(1,1))
plot(c(alldiets_mean_pupation_percentage),c(alldiets_mean_eclosion_percentage))
cor(c(alldiets_mean_pupation_percentage),c(alldiets_mean_eclosion_percentage))
hist(alldiets_mean_pupation_percentage,20)
hist(alldiets_mean_eclosion_percentage,20)

alldiets_all_percentages = rbind(cbind("Type"="Pupation",melt(alldiets_mean_pupation_percentage)),
                                 cbind("Type"="Eclosion",melt(alldiets_mean_eclosion_percentage)))
colnames(alldiets_all_percentages) <- c("Type","line","diet","percentage")

# order lines by hclust of pupation
# lineordering = hclust(dist(alldiets_mean_pupation_percentage))$order
for (orderingtype in c("none","pupation","eclosion")) {
  lineorderingpupation = names(sort(alldiets_mean_pupation_percentage[,"HPD"]))
  lineorderingeclosion = names(sort(alldiets_mean_eclosion_percentage[,"HPD"]))
  lineordering = mixedsort(rownames(alldiets_mean_eclosion_percentage))
  if (orderingtype == "none") alldiets_all_percentages$line = factor(alldiets_all_percentages$line,
                                                                     levels = lineordering)
  if (orderingtype == "pupation") alldiets_all_percentages$line = factor(alldiets_all_percentages$line,
                                                                         levels = lineorderingpupation)
  if (orderingtype == "eclosion") alldiets_all_percentages$line = factor(alldiets_all_percentages$line,
                                                                         levels = lineorderingeclosion)
  alldiets_all_percentages$diet = factor(alldiets_all_percentages$diet,
                                         levels = rev(diets))
  
  g = ggplot(alldiets_all_percentages,aes(x = line,y=diet,fill=percentage)) +
    geom_tile() + facet_wrap(~Type, ncol=1) + 
    labs(x=NULL, y=NULL, title="Percentage of pupation and eclosion by diet and line") +
    scale_fill_viridis(name="Percentage") +  theme_tufte(base_family="Helvetica") + 
    theme(axis.ticks=element_blank()) +  geom_tile(color="grey", size=0.05) +
    theme(axis.text=element_text(size=10)) + theme(panel.border=element_blank()) +
    theme(plot.title=element_text(hjust=0)) + theme(strip.text=element_text(hjust=0)) +
    theme(legend.title=element_text(size=6)) + theme(legend.title.align=1) + 
    theme(legend.text=element_text(size=6)) + theme(legend.position="bottom") + 
    theme(legend.key.size=unit(0.2, "cm")) + theme(legend.key.width=unit(1, "cm")) +
    theme(axis.text.x=element_text(size=5,angle=90,hjust=0))
  g  
  ggsave(paste0("../Figures/percentage_heatmap_ordered_by_",orderingtype,".pdf"),height=6,width=12)
  
}

# there is both a diet and line effect in both pupation and eclosion
tmp = melt(alldiets_mean_pupation_percentage)
fit = lm(value ~ X1+X2,tmp)
anova(fit)

tmp = melt(alldiets_mean_eclosion_percentage)
fit = lm(value ~ X1+X2,tmp)
anova(fit)

save(alldiets_mean_pupation_percentage,alldiets_mean_eclosion_percentage,
     alldiets_all_percentages,
     file="../Data/alldiets_mean_percentage.RData")