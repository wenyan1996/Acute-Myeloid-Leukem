library(cgdsr)
library(DT)
library(survival)
library(PAmeasures)
### 1. Data Preprocessing
## Clinical profiles
mycgds <- CGDS("http://www.cbioportal.org/")
all_TCGA_studies <- getCancerStudies(mycgds)
DT::datatable(all_TCGA_studies)
laml <- 'laml_tcga'
all_tables <- getCaseLists(mycgds, laml)
clinicaldata <- getClinicalData(mycgds, "laml_tcga_rna_seq_v2_mrna")
DT::datatable(clinicaldata,
              extensions = "FixedColumns",
              options = list(
                scrollX = TRUE,
                fixedColumns = TRUE
              ))

## Extract gene expression profiles
my_dataset <- 'laml_tcga_rna_seq_v2_mrna'
my_table <- 'laml_tcga_rna_seq_v2_mrna'
gene_set <- read.table('gene.txt')
gene <-  as.data.frame(getProfileData(mycgds, toString(gene_set[1,1]), my_dataset, my_table))
for (i in 2:nrow(gene_set)){
  g <- as.data.frame(getProfileData(mycgds, toString(gene_set[i,1]), my_dataset, my_table))
  if (nrow(g) == nrow(gene)){
    gene <- cbind(gene, g)
  }
  else {
    print(gene_set[i,1])
  }
}

## Data combination  
clinical <- clinicaldata[,c(1, 2, 4, 5, 12, 14, 15, 16, 21, 23, 24, 27, 33)]

## Remove missing value and subjects with 0 months overall survival
clinical <- merge(clinical, gene, by = "row.names")
clinical$delta <- ifelse(clinical$OS_STATUS == 'DECEASED', 1, 0)
id <- c(1:nrow(clinical))
rownames(clinical) <- id
clinical <- cbind(id, clinical)
clinical <- clinical[which(is.na(clinical$OS_MONTHS) == FALSE),]
clinical <- clinical[which(clinical$OS_MONTHS != 0),]
clinical$HISTORY_NEOADJUVANT_TRTYN <- ifelse(clinical$HISTORY_NEOADJUVANT_TRTYN == 'Yes', 1, 0)
clinical$HISTORY_OTHER_MALIGNANCY <- ifelse(clinical$HISTORY_OTHER_MALIGNANCY == 'Yes', 1, 0)
clinical$MALE <- ifelse(clinical$SEX == 'Male', 1, 0)
clinical <- clinical[which(clinical$FAB != 'Not Classified'),]

## Generate dummy variables 
dummy.matrix.generator <- function(varname){
  str.vector <- clinical[,varname]
  categories <- na.omit(unique(str.vector))
  K <- length(categories)
  N <- length(str.vector)
  
  mat <- data.frame(matrix(NA, nrow = N, ncol = K-1))
  names(mat) <- paste(varname, 1:(K-1), sep = "_")
  for (k in 1:(K-1)){
    mat[,k] <- (str.vector == categories[k]) * 1
  }
  return(mat)
}
FAB_matrix <- dummy.matrix.generator("FAB")
clinical <- cbind(clinical, FAB_matrix)

### After data preprocessing, clinical features: 17; genes: 61; 148 subjects.

### 2. Explore the data
library(ggplot2)
ggplot(clinical,aes(x = OS_MONTHS, colour = OS_STATUS)) +
  geom_density(alpha = 0.5) +
  ggtitle("Distribution of Overall Survival Time for Censored and Uncensored Data") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("Overall Survival Time (Months)") +
  scale_color_manual(labels = c("Deceased (Uncensored)", "Living (Censored)"),
                     values = c("red", "blue"),
                     name = "Status") 

table(clinical$OS_STATUS)
my.surv <- Surv(clinical$OS_MONTHS, clinical$OS_STATUS == 'DECEASED')
kmfit <- survfit(my.surv ~ 1)
summary(kmfit)
plot(kmfit)

gene <- clinical[,c(16:74)]
result <- c()
for (i in 1:ncol(gene)){
  value <- gene[,i]
  group <- ifelse(value > median(value), 'high', 'low')
  kmfit <- survfit(my.surv ~ group)
  data.survdiff <- survdiff(my.surv ~ group)
  p.value <- 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  if (p.value < 0.05){
    result[length(result) + 1] <- names(gene[i])
    print(names(gene[i]))
    print(p.value)
  }
}

## result: PIK3CA PIK3CB BCL2A1 MAP2K1 PPARD PIM1 KIT ITGAM MPO STAT5B
kmfit.2 <- survfit(Surv(OS_MONTHS, delta) ~ FAB, data = clinical)
plot(kmfit.2, col = rainbow((length(unique(clinical$FAB)))))

survdiff(Surv(OS_MONTHS, delta) ~ FAB, data = clinical)
survdiff(Surv(OS_MONTHS, delta) ~ SEX, data = clinical)
survdiff(Surv(OS_MONTHS, delta) ~ HISTORY_NEOADJUVANT_TRTYN, data = clinical)
survdiff(Surv(OS_MONTHS, delta) ~ HISTORY_OTHER_MALIGNANCY, data = clinical)

### Clinical
fit.clinical.cox <- coxph(Surv(OS_MONTHS, delta) ~ ABNORMAL_LYMPHOCYTE_PERCENT + AGE + BASOPHILS_COUNT + HISTORY_NEOADJUVANT_TRTYN + 
                      HISTORY_OTHER_MALIGNANCY + BLAST_COUNT + FRACTION_GENOME_ALTERED + MUTATION_COUNT +  PLATELET_COUNT_PRERESECTION + 
                      MALE + FAB_1 + FAB_2 + FAB_3 + FAB_4 + FAB_5 + FAB_6 + FAB_7, 
                      data = clinical, ties = 'breslow',  x = TRUE, y = TRUE)
fit.clinical.cox.best <- coxph(Surv(OS_MONTHS, delta) ~ ABNORMAL_LYMPHOCYTE_PERCENT + AGE + HISTORY_NEOADJUVANT_TRTYN + 
                            BLAST_COUNT + FRACTION_GENOME_ALTERED,
                          data = clinical, ties = 'breslow',  x = TRUE, y = TRUE)
summary(fit.clinical.cox.best)
stepAIC(fit.clinical.cox)
pam.coxph(fit.clinical.cox)  # R-squared = 0.1405, L-squared = 0.3634

fit.clinical.lognormal <- survreg(Surv(OS_MONTHS, delta) ~ ABNORMAL_LYMPHOCYTE_PERCENT + AGE + BASOPHILS_COUNT + HISTORY_NEOADJUVANT_TRTYN + 
                                    HISTORY_OTHER_MALIGNANCY + BLAST_COUNT + FRACTION_GENOME_ALTERED + MUTATION_COUNT +  PLATELET_COUNT_PRERESECTION + 
                                    MALE + FAB_1 + FAB_2 + FAB_3 + FAB_4 + FAB_5 + FAB_6 + FAB_7,  data = clinical, dist = 'lognormal', x = TRUE, y = TRUE)
summary(fit.clinical.lognormal)
stepAIC(fit.clinical.lognormal)
fit.clinical.lognormal.best <- survreg(Surv(OS_MONTHS, delta) ~ ABNORMAL_LYMPHOCYTE_PERCENT + AGE + HISTORY_NEOADJUVANT_TRTYN + 
                                     FRACTION_GENOME_ALTERED + FAB_1 + FAB_2 + FAB_3 + FAB_4 + FAB_5 + FAB_6 + FAB_7,  data = clinical, dist = 'lognormal', x = TRUE, y = TRUE)
summary(fit.clinical.lognormal.best)
pam.survreg(fit.clinical.lognormal)  # R-squared = 0.0039, L-squared = 0.9826

fit.clinical.weibull <- survreg(Surv(OS_MONTHS, delta) ~ ABNORMAL_LYMPHOCYTE_PERCENT + AGE + BASOPHILS_COUNT + HISTORY_NEOADJUVANT_TRTYN + 
                                    HISTORY_OTHER_MALIGNANCY + BLAST_COUNT + FRACTION_GENOME_ALTERED + MUTATION_COUNT +  PLATELET_COUNT_PRERESECTION + 
                                    MALE + FAB_1 + FAB_2 + FAB_3 + FAB_4 + FAB_5 + FAB_6 + FAB_7,  data = clinical, dist = 'weibull', x = TRUE, y = TRUE)
summary(fit.clinical.weibull)
stepAIC(fit.clinical.weibull)
fit.clinical.weibull.best <- survreg(Surv(OS_MONTHS, delta) ~ ABNORMAL_LYMPHOCYTE_PERCENT + AGE  + HISTORY_NEOADJUVANT_TRTYN + 
                                  BLAST_COUNT + FRACTION_GENOME_ALTERED +
                                  FAB_1 + FAB_2 + FAB_3 + FAB_4 + FAB_5 + FAB_6 + FAB_7,  data = clinical, dist = 'weibull', x = TRUE, y = TRUE)
summary(fit.clinical.weibull.best)
pam.survreg(fit.clinical.weibull)  # R-squared = 0.0031, L-squared = 0.9834

### Gene 
fit.gene.cox <- coxph(Surv(OS_MONTHS, delta) ~  PIK3CA + PIK3CB + BCL2A1 + MAP2K1 
                  + PPARD + PIM1 + KIT + ITGAM + MPO + STAT5B, data = clinical, ties = 'breslow', x = TRUE, y = TRUE)
summary(fit.gene.cox)
pam.coxph(fit.gene.cox)  # R-squared: 0.1525, L-squared = 0.4163
stepAIC(fit.gene.cox)
fit.gene.cox.best <- coxph(Surv(OS_MONTHS, delta) ~  PIK3CB + PIM1 +
                            MPO + STAT5B, data = clinical, ties = 'breslow', x = TRUE, y = TRUE)
summary(fit.gene.cox.best)
fit.gene.lognormal <- survreg(Surv(OS_MONTHS, delta) ~  PIK3CA + PIK3CB + BCL2A1 + MAP2K1 
                    + PPARD + PIM1 + KIT + ITGAM + MPO + STAT5B, data = clinical, dist = 'lognormal', x = TRUE, y = TRUE)
summary(fit.gene.lognormal)
stepAIC(fit.gene.lognormal)
fit.gene.lognormal.best <- survreg(Surv(OS_MONTHS, delta) ~  PIK3CA + PIM1 + MPO, data = clinical, dist = 'lognormal', x = TRUE, y = TRUE)
summary(fit.gene.lognormal.best)
pam.survreg(fit.gene.lognormal)   # R-squared = 0.0031, L-squared = 0.9834

fit.gene.weibull <- survreg(Surv(OS_MONTHS, delta) ~  PIK3CA + PIK3CB + BCL2A1 + MAP2K1 
                              + PPARD + PIM1 + KIT + ITGAM
                              + MPO + STAT5B, data = clinical, dist = 'weibull', x = TRUE, y = TRUE)
summary(fit.gene.weibull)
pam.survreg(fit.gene.weibull)   # R-squared = 0.0010, L-squared = 0.9855
stepAIC(fit.gene.weibull)
fit.gene.weibull.best <- survreg(Surv(OS_MONTHS, delta) ~  PIK3CB + PIM1 + MPO + STAT5B, data = clinical, dist = 'weibull',
                                 x = TRUE, y = TRUE)
summary(fit.gene.weibull.best)


### Combination 
fit.combine.cox <-  coxph(Surv(OS_MONTHS, delta) ~ ABNORMAL_LYMPHOCYTE_PERCENT + AGE + BASOPHILS_COUNT + HISTORY_NEOADJUVANT_TRTYN + 
                        HISTORY_OTHER_MALIGNANCY + BLAST_COUNT + FRACTION_GENOME_ALTERED + MUTATION_COUNT +  PLATELET_COUNT_PRERESECTION + 
                        MALE +  FAB_1 + FAB_2 + FAB_3 + FAB_4 + FAB_5 + FAB_6 + FAB_7 + PIK3CA + PIK3CB + BCL2A1 + MAP2K1 + 
                        PPARD + PIM1 + KIT  + ITGAM + MPO + STAT5B, data = clinical, ties = 'breslow', x = TRUE, y = TRUE)
summary(fit.combine.cox)
stepAIC(fit.combine.cox)
fit.combine.cox.best <-  coxph(Surv(OS_MONTHS, delta) ~ ABNORMAL_LYMPHOCYTE_PERCENT + AGE + HISTORY_NEOADJUVANT_TRTYN + 
                            FRACTION_GENOME_ALTERED + PIK3CA + PIM1 + MPO, data = clinical, ties = 'breslow', x = TRUE, y = TRUE)
summary(fit.combine.cox.best)
pam.coxph(fit.combine.cox)  # R-squared = 0.2049, L-squared = 0.3765
fit.combine.lognormal <- survreg(Surv(OS_MONTHS, delta) ~ ABNORMAL_LYMPHOCYTE_PERCENT + AGE + BASOPHILS_COUNT + HISTORY_NEOADJUVANT_TRTYN + 
                                   HISTORY_OTHER_MALIGNANCY + BLAST_COUNT + FRACTION_GENOME_ALTERED + MUTATION_COUNT +  PLATELET_COUNT_PRERESECTION + 
                                   MALE + FAB_1 + FAB_2 + FAB_3 + FAB_4 + FAB_5 + FAB_6 + FAB_7 + PIK3CA + PIK3CB + BCL2A1 + MAP2K1 
                                 + PPARD + PIM1 + KIT  + ITGAM + MPO + STAT5B, data = clinical, dist = 'lognormal', x = TRUE, y = TRUE)

stepAIC(fit.combine.lognormal)
fit.combine.lognormal.best <- survreg(Surv(OS_MONTHS, delta) ~ ABNORMAL_LYMPHOCYTE_PERCENT + AGE + FRACTION_GENOME_ALTERED + PIK3CA + MAP2K1 
                                 + PPARD + PIM1 + MPO, data = clinical, dist = 'lognormal', x = TRUE, y = TRUE)
summary(fit.combine.lognormal.best)
pam.survreg(fit.combine.lognormal) # R-squared = 0.0220, L-squared = 0.9647
fit.combine.weibull <- survreg(Surv(OS_MONTHS, delta) ~ ABNORMAL_LYMPHOCYTE_PERCENT + AGE + BASOPHILS_COUNT + HISTORY_NEOADJUVANT_TRTYN + 
                                   HISTORY_OTHER_MALIGNANCY + BLAST_COUNT + FRACTION_GENOME_ALTERED + MUTATION_COUNT +  PLATELET_COUNT_PRERESECTION + 
                                   MALE + FAB_1 + FAB_2 + FAB_3 + FAB_4 + FAB_5 + FAB_6 + FAB_7 + PIK3CA + PIK3CB + BCL2A1 + MAP2K1 
                                 + PPARD + PIM1 + KIT  + ITGAM + MPO + STAT5B, data = clinical, dist = 'weibull', x = TRUE, y = TRUE)
pam.survreg(fit.combine.weibull)  # R-squared = 0.0074, L-squared = 0.9791
stepAIC(fit.combine.weibull)
fit.combine.weibull.best <- survreg(Surv(OS_MONTHS, delta) ~ ABNORMAL_LYMPHOCYTE_PERCENT + AGE + FRACTION_GENOME_ALTERED + HISTORY_NEOADJUVANT_TRTYN
                                      + PIK3CB + PIM1 + MPO + STAT5B, data = clinical, dist = 'weibull', x = TRUE, y = TRUE)
summary(fit.combine.weibull.best)
library(MASS)
stepAIC(fit.combine.cox, direction = "both")
stepAIC(fit.combine.weibull, direction = "both")

fit.final <- coxph(Surv(OS_MONTHS, delta) ~ ABNORMAL_LYMPHOCYTE_PERCENT + AGE + HISTORY_NEOADJUVANT_TRTYN +  FRACTION_GENOME_ALTERED
                     + PIK3CA + PIM1 + MPO, data = clinical, ties = 'breslow', x = TRUE, y = TRUE)
pam.coxph(fit.final)

fit.final.lognormal <- survreg(Surv(OS_MONTHS, delta) ~ ABNORMAL_LYMPHOCYTE_PERCENT + AGE + HISTORY_NEOADJUVANT_TRTYN +  FRACTION_GENOME_ALTERED
                                + PIK3CA + PIM1 + MPO, data = clinical, dist = 'lognormal', x = TRUE, y = TRUE)
pam.survreg(fit.final.lognormal)  # R-squared = 0.0537  # L-squared = 0.9335

fit.final.weibull <- survreg(Surv(OS_MONTHS, delta) ~ ABNORMAL_LYMPHOCYTE_PERCENT + AGE + HISTORY_NEOADJUVANT_TRTYN +  FRACTION_GENOME_ALTERED
                               + PIK3CA + PIM1 + MPO, data = clinical, dist = 'weibull', x = TRUE, y = TRUE)
pam.survreg(fit.final.weibull)  # R-squared = 0.0694  # L-squared = 0.9180

library("survminer")
my.surv <- Surv(clinical$OS_MONTHS, clinical$delta)
kmfit.1 <- survfit(my.surv ~ 1)
kmfit.2 <- survfit(my.surv ~ FAB, data = clinical)
kmfit.3 <- survfit(my.surv ~ HISTORY_NEOADJUVANT_TRTYN, data = clinical)
ggsurvplot(kmfit.3, conf.int = FALSE, pval = TRUE, risk.table = TRUE, plot = TRUE, data = clinical)

pam.coxph(fit.clinical.cox.best)
pam.coxph(fit.gene.cox.best)
pam.coxph(fit.combine.cox.best)

pam.survreg(fit.clinical.lognormal.best)
pam.survreg(fit.gene.lognormal.best)
pam.survreg(fit.combine.lognormal.best)

pam.survreg(fit.clinical.weibull.best)
pam.survreg(fit.gene.weibull.best)
pam.survreg(fit.combine.weibull.best)

mg.residual <- resid(fit.combine.cox.best, type = "martingale")
cs.residual <- clinical$delta - mg.residual
fit.cs <- survfit(Surv(cs.residual, clinical$delta) ~ 1) 
H.cs   <- cumsum(fit.cs$n.event/fit.cs$n.risk)
plot(fit.cs$time, H.cs, type='s', col='blue', main = 'Cox-Snell Residual Plot', xlab = 'Residual', ylab = 'Nelson-Aalen Cum. Hazard')
abline(0, 1, col='red',  lty=2)

mg.residual <- resid(fit.combine.cox, type = "martingale")
cs.residual <- clinical$delta - mg.residual
fit.cs <- survfit(Surv(cs.residual, clinical$delta) ~ 1) 
H.cs   <- cumsum(fit.cs$n.event/fit.cs$n.risk)
plot(fit.cs$time, H.cs, type='s', col='blue', main = 'Cox-Snell Residual Plot', xlab = 'Residual', ylab = 'Nelson-Aalen Cum. Hazard')
abline(0, 1, col='red',  lty=2)
