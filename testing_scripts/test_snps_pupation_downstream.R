# Compare the p-values from the two competing models

library(scattermore)
library(ggplot2)
library(patchwork)

fig.h = 8

out = readRDS("../verona/manova_pvals_coefs.Rds")

# scatterplot of -logP values between model types
# and correlation coefficients

df = data.frame(full = -log10(out$manova_pvalPupation_batch),
                old = -log10(out$manova_pvalPupation_wolb))
corval = cor(df, use = "p", method = "spearman")[1,2]

g = ggplot(df, aes(x = old, y = full)) + 
  geom_scattermore() +
  xlab("-log(P-value) previous model") +
  ylab("-log(P-value) full model") + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_classic() + 
  ggtitle(paste0("Pupation - Spearman correlation = ", signif(corval, 3))) +
  NULL
ggsave(g, file = "../Figures/modelComparison/modelComparison_pupation.pdf",
       height = fig.h, width = fig.h)


df = data.frame(full = -log10(out$manova_pvalEclosion_batch),
                old = -log10(out$manova_pvalEclosion_wolb))
corval = cor(df, use = "p", method = "spearman")[1,2]

g = ggplot(df, aes(x = old, y = full)) + 
  geom_scattermore() +
  xlab("-log(P-value) previous model") +
  ylab("-log(P-value) full model") + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_classic() + 
  ggtitle(paste0("Eclosion - Spearman correlation = ", signif(corval, 3))) +
  NULL
ggsave(g, file = "../Figures/modelComparison/modelComparison_eclosion.pdf",
       height = fig.h, width = fig.h)

diets = colnames(out$manova_coefPupation_batch)

for (diet in diets) {
  df = data.frame(full = out$manova_coefPupation_batch[,diet],
                  old = out$manova_coefPupation_wolb[,diet])
  corval = cor(df, use = "p", method = "spearman")[1,2]
  
  g = ggplot(df, aes(x = old, y = full)) + 
    geom_scattermore() +
    xlab("Fitted coefficient previous model") +
    ylab("Fitted coefficient full model") + 
    geom_abline(slope = 1, intercept = 0) + 
    theme_classic() + 
    ggtitle(paste0("Pupation - ", diet, " - Spearman correlation = ", signif(corval, 3))) +
    NULL
  ggsave(g, file = paste0("../Figures/modelComparison/modelComparison_pupation_coef_", diet, ".pdf"),
         height = fig.h, width = fig.h)
  
  df = data.frame(full = out$manova_coefEclosion_batch[,diet],
                  old = out$manova_coefEclosion_wolb[,diet])
  corval = cor(df, use = "p", method = "spearman")[1,2]
  
  g = ggplot(df, aes(x = old, y = full)) + 
    geom_scattermore() +
    xlab("Fitted coefficient previous model") +
    ylab("Fitted coefficient full model") + 
    geom_abline(slope = 1, intercept = 0) + 
    theme_classic() + 
    ggtitle(paste0("Eclosion - ", diet, " - Spearman correlation = ", signif(corval, 3))) +
    NULL
  ggsave(g, file = paste0("../Figures/modelComparison/modelComparison_eclosion_coef_", diet, ".pdf"),
         height = fig.h, width = fig.h)
}


# repeat for HPD normalised data

out = readRDS("../verona/manova_pvals_coefs_HPD.Rds")

# scatterplot of -logP values between model types
# and correlation coefficients

df = data.frame(full = -log10(out$manova_pvalPupation_batch),
                old = -log10(out$manova_pvalPupation_wolb))
corval = cor(df, use = "p", method = "spearman")[1,2]

g = ggplot(df, aes(x = old, y = full)) + 
  geom_scattermore() +
  xlab("-log(P-value) previous model") +
  ylab("-log(P-value) full model") + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_classic() + 
  ggtitle(paste0("Pupation (HPD normalised) - Spearman correlation = ", signif(corval, 3))) +
  NULL
ggsave(g, file = "../Figures/modelComparison/modelComparison_pupationHPD.pdf",
       height = fig.h, width = fig.h)


df = data.frame(full = -log10(out$manova_pvalEclosion_batch),
                old = -log10(out$manova_pvalEclosion_wolb))
corval = cor(df, use = "p", method = "spearman")[1,2]

g = ggplot(df, aes(x = old, y = full)) + 
  geom_scattermore() +
  xlab("-log(P-value) previous model") +
  ylab("-log(P-value) full model") + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_classic() + 
  ggtitle(paste0("Eclosion (HPD normalised) - Spearman correlation = ", signif(corval, 3))) +
  NULL
ggsave(g, file = "../Figures/modelComparison/modelComparison_eclosionHPD.pdf",
       height = fig.h, width = fig.h)

diets = colnames(out$manova_coefPupation_batch)

for (diet in diets) {
  df = data.frame(full = out$manova_coefPupation_batch[,diet],
                  old = out$manova_coefPupation_wolb[,diet])
  corval = cor(df, use = "p", method = "spearman")[1,2]
  
  g = ggplot(df, aes(x = old, y = full)) + 
    geom_scattermore() +
    xlab("Fitted coefficient previous model") +
    ylab("Fitted coefficient full model") + 
    geom_abline(slope = 1, intercept = 0) + 
    theme_classic() + 
    ggtitle(paste0("Pupation (HPD normalised)  - ", diet, " - Spearman correlation = ", signif(corval, 3))) +
    NULL
  ggsave(g, file = paste0("../Figures/modelComparison/modelComparison_pupationHPD_coef_", diet, ".pdf"),
         height = fig.h, width = fig.h)
  
  df = data.frame(full = out$manova_coefEclosion_batch[,diet],
                  old = out$manova_coefEclosion_wolb[,diet])
  corval = cor(df, use = "p", method = "spearman")[1,2]
  
  g = ggplot(df, aes(x = old, y = full)) + 
    geom_scattermore() +
    xlab("Fitted coefficient previous model") +
    ylab("Fitted coefficient full model") + 
    geom_abline(slope = 1, intercept = 0) + 
    theme_classic() + 
    ggtitle(paste0("Eclosion (HPD normalised)  - ", diet, " - Spearman correlation = ", signif(corval, 3))) +
    NULL
  ggsave(g, file = paste0("../Figures/modelComparison/modelComparison_eclosionHPD_coef_", diet, ".pdf"),
         height = fig.h, width = fig.h)
}

# example barplots for the groups

df_example_A = data.frame(
  phenotype = c(1, 0.75, 0.75, 0.5),
  snp = factor(c(0,1,0,1)),
  diet = factor(c("HPD", "HPD", "HSD", "HSD"))
)

gA = ggplot(df_example_A, aes(x = snp, y = phenotype)) + 
  geom_point(shape = "_", size = 30) + 
  facet_grid(~diet) +
  ylim(c(0,1)) + 
  xlab("SNP") +
  ylab("Phenotype value") +
  ggtitle("Diet independent") +
  theme_classic()


df_example_B = data.frame(
  phenotype = c(0, 0, -0.25, -0.25),
  snp = factor(c(0,1,0,1)),
  diet = factor(c("HPD", "HPD", "HSD", "HSD"))
)

gB = ggplot(df_example_B, aes(x = snp, y = phenotype)) + 
  geom_point(shape = "_", size = 30) + 
  facet_grid(~diet) +
  ylim(c(-0.5,0.5)) + 
  xlab("SNP") +
  ylab("Normalised phenotype value") +
  ggtitle(" ") +
  theme_classic()


df_example_C = data.frame(
  phenotype = c(1, 0.75, 0.8, 0.8),
  snp = factor(c(0,1,0,1)),
  diet = factor(c("HPD", "HPD", "HSD", "HSD"))
)

gC = ggplot(df_example_C, aes(x = snp, y = phenotype)) + 
  geom_point(shape = "_", size = 30) + 
  facet_grid(~diet) +
  ylim(c(0,1)) + 
  xlab("SNP") +
  ylab("Phenotype value") +
  ggtitle("Diet dependent") +
  theme_classic()


df_example_D = data.frame(
  phenotype = c(0, 0, -0.2, 0.05),
  snp = factor(c(0,1,0,1)),
  diet = factor(c("HPD", "HPD", "HSD", "HSD"))
)

gD = ggplot(df_example_D, aes(x = snp, y = phenotype)) + 
  geom_point(shape = "_", size = 30) + 
  facet_grid(~diet) +
  ylim(c(-0.5,0.5)) + 
  xlab("SNP") +
  ylab("Normalised phenotype value") +
  ggtitle(" ") +
  theme_classic()

g = gA + gB + gC + gD + plot_layout(nrow = 2, ncol = 2, byrow = TRUE)

g
ggsave(g, file = "../Figures/modelComparison/example_HPD_normalisation.pdf",
       height = 6, width = 8)

numList = readRDS("../verona/numList.Rds")
pdf("../Figures/modelComparison/DGRP_numberAlleles.pdf",
    height = 8, width = 5)
par(mfrow = c(2,1))
hist(numList[[2]], 100, main = "Number of alternate alleles in the DGRP",
     xlab = "", ylab = "Proportion", prob = TRUE)
abline(v = 5, lty = 2)
hist(numList[[2]][numList[[2]] > 5], 100,
     main = "Number of alternate alleles in the DGRP\n(greater than 5)",
     xlab = "", ylab = "Proportion", prob = TRUE, xlim = c(5, 196))
dev.off()

