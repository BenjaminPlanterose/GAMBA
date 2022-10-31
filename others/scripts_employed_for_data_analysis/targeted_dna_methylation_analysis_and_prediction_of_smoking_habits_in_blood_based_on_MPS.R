############################################################################
############################################################################
###########                                                      ###########
###########             Author: Benjamin Planterose              ###########
###########                                                      ###########
###########        Erasmus MC University Medical Centre          ###########
###########               Rotterdam, The Netherlands             ###########
###########                                                      ###########
###########                                                      ###########
############################################################################
############################################################################

## Load libraries
library(data.table)
library(gplots)
library(ggplot2)
library(ggridges)
library(ggtern)
library(RColorBrewer)
library(scales)
library(corrplot)
library(MASS)
library(nnet)
library(caret)
library(minpack.lm)
library(ROCR)
library(HandTill2001)
library(MLmetrics)
library(effectsize)
library(UMtools)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(UMtools)

## Load functions
process_fread <- function(mat)
{
  CpGs = mat$rn
  mat = as.matrix(mat[,-1])
  rownames(mat) = CpGs
  return(mat)
}

process.names <- function(X, run)
{
  Nam = colnames(X)
  Nam = paste(run, sapply(strsplit(Nam, split = "_"), function(x) x[3]), sep = "")
  colnames(X) = Nam
  return(X)
}

inverse.model <- function(k1, k2, a, b, val)
{
  f <- function(x)
  {
    (a*exp(-k1*x) + b*exp(-k2*x) - val)
  }
  uniroot(f = f, interval = c(0,1))$root
}

f1 <- function(data, lev = NULL, model = NULL) 
{
  confusionMatrix(data = table(data$pred, data$obs), mode = 'everything')$byClass['F1']
}

macrof1 <- function(data, lev = NULL, model = NULL) 
{
  c(F1 = F1_Score_macro_weighted(data$obs, data$pred))
}

####################### General #######################

breaks=seq(0, 1, 0.05)
my_palette <- colorRampPalette(c("ivory", "cyan4"))(n = length(breaks) - 1)
LEV = c("cg05575921", "cg13039251", "cg03636183", "cg12803068", "cg22132788",
        "cg06126421", "cg21566642", "cg23576855", "cg15693572", "cg05951221",
        "cg01940273", "cg12876356", "cg09935388")

setwd("<not_shown>")
cg_pos = fread("cg_position.txt")
cg_pos = cg_pos[order(cg_pos$cg),]
POS_450K = paste(cg_pos$cg, as.numeric(cg_pos$position)+2, as.numeric(cg_pos$position)+2, sep = "_") # +2 because of NN, beginning amplicon

setwd("<not_shown>")
load(file = "Two_category_model.rda")
M2cat
load(file = "Three_category_model.rda")
M3cat

####################### Comparing runs #######################

# Run 1
setwd("<not_shown>")
counts = fread("counts.txt", header = F)$V1
Samp = counts[seq(1, length(counts), 2)]
counts = as.numeric(counts[seq(2, length(counts), 2)])
names(counts) = Samp
barplot(log10(counts+1)); abline(h = 3, lty = 2, col = 'red2')
(failed = Samp[log10(counts+1) < 3]) # S45, S58, S73, S75
counts = counts[log10(counts+1) > 3] # remove failed barcounts
sum(counts)/10^6 # 30.29352 total reads
mean(counts) # 172122.3 read pairs per sample
mean(counts)/13 # 13240.18 read pairs per marker per sample
setwd("<not_shown>")
METH_cov = process_fread(fread("2022-01-24_coverage.txt", header = T))
sum(POS_450K %in% rownames(METH_cov)) # 13
failed = unique(sapply(strsplit(failed, split = "_"), function(x) x[2]))
remove = unlist(sapply(1:length(failed), function(x) grep(paste("_", failed[x], "_", sep = ""), colnames(METH_cov))))
cov.m1 = melt(METH_cov[POS_450K, !(1:ncol(METH_cov) %in% remove)])
cov.m1$Var1 = sapply(strsplit(as.vector(cov.m1$Var1), split = '_'), function(x) x[1])
ggplot(data = cov.m1, mapping = aes(x = Var1, y = log10(value+1), fill = Var1)) + geom_violin() +
labs(title = "Run 1", x = "CpG", y = "log10(read_count + 1)") + ylim(0, 6) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        plot.title = element_text(hjust = 0.5))

# Run 2
setwd("<not_shown>")
counts = fread("counts.txt", header = F)$V1
Samp = counts[seq(1, length(counts), 2)]
counts = as.numeric(counts[seq(2, length(counts), 2)])
names(counts) = Samp
barplot(log10(counts+1)); abline(h = 3, lty = 2, col = 'red2')
(failed = Samp[log10(counts+1) < 3]) # S19, S51, S54, S62, S64, S78
counts = counts[log10(counts+1) > 3] # remove failed barcounts
sum(counts)/10^6 # 29.42786 M total reads
mean(counts) # 193604.4 read pairs per sample
mean(counts)/13 # 14892.64 read pairs per marker per sample
setwd("<not_shown>")
METH_cov = process_fread(fread("2022-01-25_coverage.txt", header = T))
sum(POS_450K %in% rownames(METH_cov)) # 13
failed = unique(sapply(strsplit(failed, split = "_"), function(x) x[2]))
remove = unlist(sapply(1:length(failed), function(x) grep(paste("_", failed[x], "_", sep = ""), colnames(METH_cov))))
cov.m2 = melt(METH_cov[POS_450K, !(1:ncol(METH_cov) %in% remove)])
cov.m2$Var1 = sapply(strsplit(as.vector(cov.m2$Var1), split = '_'), function(x) x[1])

p2 = ggplot(data = cov.m2, mapping = aes(x = Var1, y = log10(value+1), fill = Var1)) + geom_violin() +
  labs(title = "Run 2", x = "CpG", y = "log10(read_count + 1)") + ylim(0, 6) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        plot.title = element_text(hjust = 0.5))

# Run 3
setwd("<not_shown>")
counts = fread("counts.txt", header = F)$V1
Samp = counts[seq(1, length(counts), 2)]
counts = as.numeric(counts[seq(2, length(counts), 2)])
names(counts) = Samp
barplot(log10(counts+1)); abline(h = 3, lty = 2, col = 'red2')
(failed = Samp[log10(counts+1) < 4.6]) # S1, S2
counts = counts[log10(counts+1) > 4.6] # remove failed barcounts
sum(counts)/10^6 # 37.7425 M total reads
mean(counts) # 200758 read pairs per sample
mean(counts)/13 # 15442.92 read pairs per marker per sample
setwd("<not_shown>")
METH_cov = process_fread(fread("2022-01-26_coverage.txt", header = T))
sum(POS_450K %in% rownames(METH_cov)) # 13
failed = unique(sapply(strsplit(failed, split = "_"), function(x) x[2]))
remove = unlist(sapply(1:length(failed), function(x) grep(paste("_", failed[x], "_", sep = ""), colnames(METH_cov))))
cov.m3 = melt(METH_cov[POS_450K, !(1:ncol(METH_cov) %in% remove)])
cov.m3$Var1 = sapply(strsplit(as.vector(cov.m3$Var1), split = '_'), function(x) x[1])

p3 = ggplot(data = cov.m3, mapping = aes(x = Var1, y = log10(value+1), fill = Var1)) + geom_violin() +
  labs(title = "Run 3", x = "CpG", y = "log10(read_count + 1)") + ylim(0, 6) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        plot.title = element_text(hjust = 0.5))

colnames(cov.m1) = colnames(cov.m2) = colnames(cov.m3) = c('CpG', 'file', 'count')
cov.m1$run = 'run1'; cov.m2$run = 'run2'; cov.m3$run = 'run3'
cov.m = Reduce(rbind, list(cov.m1, cov.m2, cov.m3))

p = ggplot(data = cov.m, aes(x = log10(count+1), y = run, fill = CpG)) + geom_density_ridges() + facet_wrap(~CpG, ncol = 3) + 
  theme_minimal() + xlab('') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position= 'none') + 
  scale_fill_manual(values = c(brewer.pal(8, 'Dark2'), brewer.pal(5, 'Accent')))
setwd("<not_shown>")
ggsave(filename = 'C_count.tiff', plot = p, device = 'tiff', width = 6.7, height = 6.2, units = 'in', dpi = 300)

####################### pooling METH #######################

# Read METH
setwd("<not_shown>")
METH1 = process_fread(fread("2022-01-24_meth.txt", header = T))/100
setwd("<not_shown>")
METH2 = process_fread(fread("2022-01-25_meth.txt", header = T))/100
setwd("<not_shown>")
METH3 = process_fread(fread("2022-01-26_meth.txt", header = T))/100

# process names and filter for 450K pos
METH1 = process.names(METH1[POS_450K,], "")
METH2 = process.names(METH2[POS_450K,], "")
METH3 = process.names(METH3[POS_450K,], "")

# Remove failed barcodes
setwd("<not_shown>")
counts = fread("counts.txt", header = F)$V1
Samp = counts[seq(1, length(counts), 2)]
counts = as.numeric(counts[seq(2, length(counts), 2)]); names(counts) = Samp
failed = Samp[log10(counts+1) < 3] # S45, S58, S73, S75
(remove = unique(sapply(strsplit(failed, split = '_'), function(x) x[2])))
METH1 = METH1[, !(colnames(METH1) %in% remove)]
#
setwd("<not_shown>")
counts = fread("counts.txt", header = F)$V1
Samp = counts[seq(1, length(counts), 2)]
counts = as.numeric(counts[seq(2, length(counts), 2)]); names(counts) = Samp
failed = Samp[log10(counts+1) < 3] # S19, S51, S54, S62, S64, S78
(remove = unique(sapply(strsplit(failed, split = '_'), function(x) x[2])))
METH2 = METH2[, !(colnames(METH2) %in% remove)]
#
setwd("<not_shown>")
counts = fread("counts.txt", header = F)$V1
Samp = counts[seq(1, length(counts), 2)]
counts = as.numeric(counts[seq(2, length(counts), 2)]); names(counts) = Samp
failed = Samp[log10(counts+1) < 4.6] # character(0)
(remove = unique(sapply(strsplit(failed, split = '_'), function(x) x[2])))
METH3 = METH3[, !(colnames(METH3) %in% remove)]

# Read phenotype
setwd("<not_shown>")
phenotype = as.data.frame(fread("ERF_phenotypes.txt"))
phenotype$sequencing_index2 = paste('S', as.integer(sapply(strsplit(phenotype$sequencing_index, split = 'UDI'), function(x) x[2])), sep = '')

# Remove STDs
index = phenotype[(phenotype$sequencing_run == 1) & startsWith(phenotype$ERF_ID, 'STD'),]$sequencing_index2
index = index[index %in% colnames(METH1)]
STD1 = METH1[, index]
colnames(STD1) = phenotype[(phenotype$sequencing_run == 1) & (phenotype$sequencing_index2 %in% index),]$ERF_ID
METH1 = METH1[, !(colnames(METH1) %in% index)]
index = phenotype[(phenotype$sequencing_run == 2) & startsWith(phenotype$ERF_ID, 'STD'),]$sequencing_index2
index = index[index %in% colnames(METH2)]
STD2 = METH2[, index]
colnames(STD2) = phenotype[(phenotype$sequencing_run == 2) & (phenotype$sequencing_index2 %in% index),]$ERF_ID
METH2 = METH2[, !(colnames(METH2) %in% index)]

# Recode METH with ERF_IDs
index = match(colnames(METH1), phenotype[phenotype$sequencing_run == 1,]$sequencing_index2)
sum(is.na(index)) # 0
colnames(METH1) = phenotype[phenotype$sequencing_run == 1,][index,]$ERF_ID
index = match(colnames(METH2), phenotype[phenotype$sequencing_run == 2,]$sequencing_index2)
sum(is.na(index)) # 0
colnames(METH2) = phenotype[phenotype$sequencing_run == 2,][index,]$ERF_ID
index = match(colnames(METH3), phenotype[phenotype$sequencing_run == 3,]$sequencing_index2)
sum(is.na(index)) # 0
colnames(METH3) = phenotype[phenotype$sequencing_run == 3,][index,]$ERF_ID
METH = Reduce(cbind, list(METH1, METH2, METH3))
rownames(METH) = sapply(strsplit(rownames(METH), split = "_"), function(x) x[1])
which(table(colnames(METH)) > 1) # 2931 3083 3869  417 4552  653 

####################### visualizing STD #######################

rownames(STD1) = sapply(strsplit(rownames(STD1), split = "_"), function(x) x[1])
rownames(STD2) = sapply(strsplit(rownames(STD2), split = "_"), function(x) x[1])
setwd("<not_shown>")
UMtools::export_bigmat(bigmat = round(rbind(t(STD1), t(STD2)), 5), 'STDs.txt')

STD = rbind(cbind(melt(t(STD1)), std = 'run 1'), cbind(melt(t(STD2)), std = 'run 2'))
head(STD)
STD$ratio = as.numeric(sapply(strsplit(as.vector(STD$Var1), split = ' '), function(x) x[2]))
head(STD)
ggplot(data = STD, aes(x = ratio, y = value, col = Var2)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_smooth(method = "loess", span = 1) + 
  theme_minimal() + geom_jitter() +  facet_wrap(~Var2) + ylim(-0.15,1.1) + 
  labs(title = "Standards") +
  xlab("Expected Methylation Ratio") + ylab("Observed Methylation Ratio") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none",
        plot.title = element_text(hjust = 0.5))

pooling_info = matrix(c('cg05575921', 2, 'cg06126421', 3, 'cg22132788', 3, 
                        'cg12803068', 3, 'cg13039251', 2, 'cg01940273', 1, 
                        'cg03636183', 1, 'cg23576855', 2, 'cg12876356', 1,
                        'cg15693572', 2, 'cg05951221', 3, 'cg09935388', 3, 'cg21566642', 1), ncol = 2, byrow = T)
pooling_info = as.data.frame(pooling_info)
pool1 = pooling_info$V1[pooling_info$V2 == '1']
pool2 = pooling_info$V1[pooling_info$V2 == '2']
pool3 = pooling_info$V1[pooling_info$V2 == '3']
STD$pool = NA
STD$pool[STD$Var2 %in% pool1] = 'pool_1'
STD$pool[STD$Var2 %in% pool2] = 'pool_2'
STD$pool[STD$Var2 %in% pool3] = 'pool_3'

p = ggplot(data = STD, aes(x = ratio, y = value, col = std)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5) +
  geom_smooth(method = "loess", span = 1, alpha = 0.5) + 
  theme_minimal() + geom_jitter(alpha = 0.5) +  facet_wrap(pool~Var2, drop = T) + 
  labs(title = "Standards") +
  xlab("Expected Methylation Ratio") + ylab("Observed Methylation Ratio") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5)) + scale_y_continuous(breaks = seq(0, 1, by = 0.5))
ggsave(filename = 'STDs.tiff', plot = p, device = 'tiff', width = 6.7, height = 6.2, units = 'in', dpi = 300)

cor = sapply(1:length(LEV), function(x) cor(STD$ratio[STD$Var2 == LEV[x]], STD$value[STD$Var2 == LEV[x]], method = 'spearman'))
names(cor) = LEV
# cg05575921 cg13039251 cg03636183 cg12803068 cg22132788 cg06126421 cg21566642 cg23576855 cg15693572 cg05951221 cg01940273 
# 0.9901598  0.9950860  0.9950860  0.9913913  0.9507505  0.9950860  0.9803075  0.9803075  0.9950860  0.9618343  0.9950860 
# cg12876356 cg09935388 
# 0.9174988  0.9852336 

####################### Correcting amplification bias #######################

############### Train amplification bias model

# Parameter estimation
std_data = cbind(STD1, STD2)
x = sapply(strsplit(colnames(std_data), split = '[ ]'), function(x) as.numeric(x[2]))
max_try = 5
set.seed(1)
results = list()
for(i in 1:nrow(std_data))
{
  print(i)
  y = std_data[i,]
  mod = tryCatch(nlsLM(formula = y ~ a*exp(-k1*x) + b*exp(-k2*x),
                       start = list(k1 = rnorm(1), k2 = rnorm(1), a = rnorm(1), b = rnorm(1)), control = list(maxiter = 100)), error = function(x) NA)
  j = 1
  while(is.na(mod)[1] & j < max_try)
  {
    mod = tryCatch(nlsLM(formula = y ~ a*exp(-k1*x) + b*exp(-k2*x),
                         start = list(k1 = rnorm(1), k2 = rnorm(1), a = rnorm(1), b = rnorm(1)), control = list(maxiter = 100)), error = function(x) NA)
    j = j + 1
  }
  results[[i]] = mod
}

# Inverse problem
COEF = sapply(results, coef)
colnames(COEF) = rownames(std_data)
range = seq(0, 1, 0.01)
RES = sapply(1:nrow(std_data), function(y) sapply(range, function(x) tryCatch(inverse.model(COEF[1,y], COEF[2,y], 
                                                                                COEF[3,y], COEF[4,y], x), error = function(e) NA)))
# Visualization
m.RES = melt(RES)
m.RES$Var2 = rownames(std_data)[m.RES$Var2]
m.RES$Var1 = range[m.RES$Var1]
colnames(m.RES) = c('value', 'Var2', 'ratio')

# Visualize fit
p = ggplot(data = STD, aes(x = ratio, y = value, col = Var2)) + 
  geom_line(data = m.RES) +
  theme_minimal() + geom_jitter(alpha = 0.5) +  facet_wrap(~Var2, drop = T) + 
  labs(title = "Standards") +
  xlab("Expected Methylation Ratio") + ylab("Observed Methylation Ratio") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none",
        plot.title = element_text(hjust = 0.5)) + scale_y_continuous(breaks = seq(0, 1, by = 0.5))

setwd("<not_shown>")
ggsave(filename = 'amplification_bias.tiff', plot = p, device = 'tiff', width = 6.7, height = 6.2, units = 'in', dpi = 300)
UMtools::export_bigmat(COEF, 'coefficients.txt')

############### Correct data

# Compute corrected data
METH_corrected = sapply(1:nrow(METH), function(y) sapply(METH[y,], function(x) tryCatch(inverse.model(COEF[1,y], COEF[2,y], 
                                                                                                COEF[3,y], COEF[4,y], x), error = function(e) NA)))
METH_corrected = t(METH_corrected)
rownames(METH_corrected) = rownames(METH)

# Out of prediction range
sum(is.na(METH_corrected)) # 128
METH_corrected[is.na(METH_corrected)] = round(METH)[is.na(METH_corrected)]

# Visualization
df = cbind(melt(METH), melt(METH_corrected)$value)
colnames(df) = c('CpG', 'ID', 'uncorrected_meth', 'corrected_meth')
dfB = m.RES
colnames(dfB) = c('uncorrected_meth', 'CpG', 'corrected_meth')

p= ggplot(data = df, aes(x = uncorrected_meth, y = corrected_meth, col = CpG)) + geom_jitter(alpha = 0.5) +  
  facet_wrap(~CpG, drop = T) + geom_line(data = dfB) +
  theme_minimal() +  facet_wrap(~CpG, drop = T) +
  labs(title = "Samples") + geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.5) +
  xlab("Uncorrected Methylation Ratio") + ylab("Corrected Methylation Ratio") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none",
        plot.title = element_text(hjust = 0.5)) + scale_y_continuous(breaks = seq(0, 1, by = 0.5))

setwd("<not_shown>")
ggsave(filename = 'corrections.tiff', plot = p, device = 'tiff', width = 6.7, height = 6.2, units = 'in', dpi = 300)

####################### visualizing METH - Heatmap #######################

# Heatmaps
y = phenotype[match(colnames(METH), phenotype$ERF_ID),]$smoking_3cat
X_uncor = METH[, y != '']
X_cor = METH_corrected[, y != '']
y = y[y != '']
y[y %in% c('current_smoker', 'past_smoker')] = 'past or current smoker'; y[y == 'non_smoker'] = 'non-smoker'
col_smoke = factor(y, levels = c('non-smoker', 'past or current smoker'))
levels(col_smoke) = c('bisque1', 'dimgray'); col_smoke = as.character(col_smoke)

setwd("<not_shown>")
tiff(filename = 'betauncorrected_smoke.tiff', width = 6, height = 6, units = 'in', res = 300)
heatmap.2(t(na.omit(t(X_uncor))), trace = "n", breaks= breaks, col = my_palette, margins = c(8,8), ColSideColors = col_smoke, 
          key.title = '', key.xlab = '')
dev.off()
tiff(filename = 'betacorrected_smoke.tiff', width = 6, height = 6, units = 'in', res = 300)
heatmap.2(t(na.omit(t(X_cor))), trace = "n", breaks= breaks, col = my_palette, margins = c(8,8), ColSideColors = col_smoke, 
          key.title = '', key.xlab = '')
dev.off()

# 3 cat- corrected
y = phenotype[match(colnames(METH), phenotype$ERF_ID),]$smoking_3cat
X_uncor = METH[, y != '']
X_cor = METH_corrected[, y != '']
y = y[y != '']
y[y == 'non_smoker'] = 'Never smoker'; y[y == 'current_smoker'] = 'Current smoker'; y[y == 'past_smoker'] = 'Past smoker'
col_smoke = factor(y)
levels(col_smoke) = c('bisque1', 'dimgray', 'gray82'); col_smoke = as.character(col_smoke)

tiff(filename = 'betacorrected_smoke_3cat.tiff', width = 6, height = 6, units = 'in', res = 300)
heatmap.2(t(na.omit(t(X_cor))), trace = "n", breaks= breaks, col = my_palette, margins = c(8,8), ColSideColors = col_smoke, 
          key.title = '', key.xlab = '')
dev.off()

tiff(filename = 'betauncorrected_smoke_3cat.tiff', width = 6, height = 6, units = 'in', res = 300)
heatmap.2(t(na.omit(t(X_uncor))), trace = "n", breaks= breaks, col = my_palette, margins = c(8,8), ColSideColors = col_smoke, 
          key.title = '', key.xlab = '')
dev.off()

####################### visualizing METH - Boxplot #######################

####### UNCORRECTED

# 2 Categories
df = as.data.frame(t(METH))
y = phenotype[match(colnames(METH), phenotype$ERF_ID),]$smoking_3cat
y[y == ''] = NA; y[y %in% c('current_smoker', 'past_smoker')] = 'past or current smoker'; y[y == 'non_smoker'] = 'non-smoker'
df$smoke = factor(y)
df = melt(df, id.vars = "smoke")
df$variable = as.character(df$variable)
df$variable = factor(df$variable, levels = LEV)
p = ggplot(data = na.omit(df), aes(x = variable, y = value, fill = smoke)) + geom_boxplot() + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values = c('bisque1', 'dimgray')) + xlab('')
setwd("<not_shown>")
ggsave(filename = '2cat_uncorrected.tiff', plot = p, device = 'tiff', width = 6.7, height = 6.2, units = 'in', dpi = 300)

# 3 Categories
df = as.data.frame(t(METH))
y = phenotype[match(colnames(METH), phenotype$ERF_ID),]$smoking_3cat; y[y == ''] = NA
df$smoke = factor(y, level = c('non_smoker', 'past_smoker', 'current_smoker'))
df = melt(df, id.vars = "smoke")
df$variable = as.character(df$variable)
df$variable = factor(df$variable, levels = LEV)
p = ggplot(data = na.omit(df), aes(x = variable, y = value, fill = smoke)) + geom_boxplot() + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = c('bisque1', 'gray82', 'dimgray')) + xlab("")
setwd("<not_shown>")
ggsave(filename = '3cat_uncorrected.tiff', plot = p, device = 'tiff', width = 6.7, height = 6.2, units = 'in', dpi = 300)

####### CORRECTED

# 2 Categories
df = as.data.frame(t(METH_corrected))
y = phenotype[match(colnames(METH_corrected), phenotype$ERF_ID),]$smoking_3cat
y[y == ''] = NA; y[y %in% c('current_smoker', 'past_smoker')] = 'past or current smoker'; y[y == 'non_smoker'] = 'non-smoker'
df$smoke = factor(y)
df = melt(df, id.vars = "smoke")
df$variable = as.character(df$variable)
df$variable = factor(df$variable, levels = LEV)
p = ggplot(data = na.omit(df), aes(x = variable, y = value, fill = smoke)) + geom_boxplot() + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = c('bisque1', 'dimgray')) + xlab("")
setwd("<not_shown>")
ggsave(filename = '2cat_corrected.tiff', plot = p, device = 'tiff', width = 6.7, height = 6.2, units = 'in', dpi = 300)

# 3 Categories
df = as.data.frame(t(METH_corrected))
y = phenotype[match(colnames(METH_corrected), phenotype$ERF_ID),]$smoking_3cat; y[y == ''] = NA
df$smoke = factor(y, level = c('non_smoker', 'past_smoker', 'current_smoker'))
df = melt(df, id.vars = "smoke")
df$variable = as.character(df$variable)
df$variable = factor(df$variable, levels = LEV)
p = ggplot(data = na.omit(df), aes(x = variable, y = value, fill = smoke)) + geom_boxplot() + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = c('bisque1', 'gray82', 'dimgray')) + xlab("")
setwd("<not_shown>")
ggsave(filename = '3cat_corrected.tiff', plot = p, device = 'tiff', width = 6.7, height = 6.2, units = 'in', dpi = 300)

####################### Relationship with chronological age #######################

X = as.data.frame(t(METH))
X$age = phenotype[match(colnames(METH), phenotype$ERF_ID),]$age
X$smoke = phenotype[match(colnames(METH), phenotype$ERF_ID),]$smoking_3cat
df = melt(X, id.vars = c('age', 'smoke'))
df$smoke = factor(df$smoke, c('non_smoker', 'past_smoker', 'current_smoker'))

p = ggplot(data = na.omit(df), aes(x = age, y = value, col = smoke)) + geom_point(alpha = 0.8) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_wrap(~variable) +
  scale_color_manual(values = c('bisque1', 'gray82', 'dimgray')) + xlab("Age (y)") + geom_smooth(method = 'lm', linetype = 2, col = 'red2', alpha = 0.7)

setwd("<not_shown>")
ggsave(filename = 'age.tiff', plot = p, device = 'tiff', width = 6.7, height = 6.2, units = 'in', dpi = 300)

cor(X[,c(colnames(X)[1:13], 'age')])[,'age']
# cg01940273  cg03636183  cg05575921  cg05951221  cg06126421  cg09935388  cg12803068  cg12876356  cg13039251  cg15693572 
# -0.27062299 -0.21991727 -0.19050967 -0.24481546 -0.37216599 -0.17422216 -0.01546472 -0.14591024 -0.41049773 -0.30860228 
# cg21566642  cg22132788  cg23576855
# -0.24056138 -0.03051519 -0.18084996

pval = sapply(1:(ncol(X)-2), function(x) summary(lm(X$age ~ X[,x]))$coefficients[2,'Pr(>|t|)'])
pval = p.adjust(pval, 'bonferroni')
names(pval) = colnames(X)[1:(ncol(X)-2)]
#   cg01940273   cg03636183   cg05575921   cg05951221   cg06126421   cg09935388   cg12803068   cg12876356   cg13039251   cg15693572 
# 2.668287e-04 7.605005e-03 3.879321e-02 1.603451e-03 3.217197e-08 8.711007e-02 1.000000e+00 3.052684e-01 4.240867e-10 1.340091e-05 
# cg21566642   cg22132788   cg23576855 
# 2.116939e-03 1.000000e+00 6.317236e-02 

mod = lm(age ~ ., data = X)
summary(mod) # Adjusted R-squared:  0.3248 
plot(X$age, predict(mod), xlim = c(0, 90), ylim = c(0,90)); abline(0,1, lty = 2)
cor(X$age, predict(mod)) # 0.5818663
mean(abs(X$age- predict(mod))) # MAE = 10.61964

####################### Replicate analysis #######################

col_vec = c(brewer.pal(8, 'Dark2'), brewer.pal(5, 'Accent'))

# Point-of-view on technical replicates
rep = names(which(table(colnames(METH)) > 1))
indices = sapply(1:length(rep), function(x) which(colnames(METH) == rep[x]))
plot(0, 1, type = 'n', ylim = c(0,1), xlim = c(0,1), main = '6 technical replicates', 
     xlab= 'beta_rep1', ylab = 'beta_rep2'); grid()
points(METH[,indices[1,]], METH[,indices[2,]], cex = 2, pch = 19,
       col = alpha(rep(palette.colors(length(rep), 'Dark2'), each = nrow(METH)), 0.6)); abline(0,1, lty = 2, col = 'gray')

# Point-of-view on CpGs
setwd("<not_shown>")
tiff(filename = 'tech_rep.tiff', width = 6, height = 6, units = 'in', res = 300)
plot(0, 1, type = 'n', ylim = c(0,1), xlim = c(0,1), main = '6 technical replicates', 
     xlab= 'beta_rep1', ylab = 'beta_rep2'); grid()
points(METH[,indices[1,]], METH[,indices[2,]], cex = 2, pch = 19,
       col = alpha(rep(col_vec, by = ncol(METH)), 0.6)); abline(0,1, lty = 2, col = 'gray')
legend('topleft', legend = rownames(METH), fill = col_vec, bty = 'n')
dev.off()

cor(as.numeric(METH[,indices[1,]]), as.numeric(METH[,indices[2,]])) # 0.9834162

setwd("<not_shown>")
UMtools::export_bigmat(bigmat = round(rbind(t(METH[,indices[1,]]), t(METH[,indices[2,]])), 5), 'REP_uncorrected.txt')

setwd("<not_shown>")
UMtools::export_bigmat(bigmat = round(rbind(t(METH_corrected[,indices[1,]]), t(METH_corrected[,indices[2,]])), 5), 'REP_corrected.txt')

######## CORRECTED

# Point-of-view on CpGs
setwd("<not_shown>")
tiff(filename = 'tech_rep_CORRECTED.tiff', width = 6, height = 6, units = 'in', res = 300)
rep = names(which(table(colnames(METH_corrected)) > 1))
indices = sapply(1:length(rep), function(x) which(colnames(METH_corrected) == rep[x]))
plot(0, 1, type = 'n', ylim = c(0,1), xlim = c(0,1), main = '6 technical replicates', 
     xlab= 'beta_rep1', ylab = 'beta_rep2'); grid()
points(METH_corrected[,indices[1,]], METH_corrected[,indices[2,]], cex = 2, pch = 19,
       col = alpha(rep(col_vec, by = ncol(METH)), 0.6)); abline(0,1, lty = 2, col = 'gray')
legend('topleft', legend = rownames(METH), fill = col_vec, bty = 'n')
dev.off()

####################### mQTL #######################

plot(density(METH['cg23576855',]))
smoke = phenotype[match(colnames(METH), phenotype$ERF_ID),]$smoking_3cat
smoke = factor(as.numeric(smoke == 'current_smoker'))
levels(smoke) = c('bisque1', 'dimgray'); smoke = as.character(smoke)

plot(METH['cg05575921',], METH['cg23576855',], xlim = c(0,1), ylim = c(0,1), col = alpha(smoke, 0.6), pch = 19)

data(annot_450K)
lapply(annot_450K, function(x) x[x$cg == 'cg23576855',])
# $CpG_SNP
# SNP_chr SNP_pos        rs REF ALT QUAL FILTER
# 1:       5  373300 rs6869832   G   A    .      .
# INFO
# 1: RS=6869832;RSPOS=373300;dbSNPBuildID=116;SSR=0;SAO=0;VP=0x050128080005110036000100;GENEINFO=AHRR:57491;WGT=1;VC=SNV;PM;PMC;SLO;INT;ASP;G5;KGPhase1;KGPhase3;CAF=0.9389,0.0611;COMMON=1;TOPMED=0.91835499490316004,0.08164500509683995
# CpG_chr CpG_start CpG_end         cg
# 1:       5    373298  373300 cg23576855

rep = names(which(table(colnames(METH)) > 1))
indices = sapply(1:length(rep), function(x) which(colnames(METH) == rep[x]))
plot(METH['cg23576855',indices[1,]], METH['cg23576855',indices[2,]], ylim = c(0,1), xlim = c(0,1), pch = 19,
     xlab= 'beta_rep1', ylab = 'beta_rep2'); abline(0,1, lty = 2, col = 'gray')

####################### Co-methylation between fragments #######################

# AHRR - cg05575921, cg23576855
smoke = phenotype
plot(METH['cg05575921',], METH['cg23576855',], xlim = c(0,1), ylim = c(0,1), pch = 19, col = alpha('black', 0.5))

# GFI1 - cg12876356, cg09935388
plot(METH['cg12876356',], METH['cg09935388',], xlim = c(0,1), ylim = c(0,1), pch = 19, col = alpha('black', 0.5))

# MYO1G - cg22132788, cg12803068
plot(METH['cg22132788',], METH['cg12803068',], xlim = c(0,1), ylim = c(0,1), pch = 19, col = alpha('black', 0.5))

# IG-chr2 - cg01940273, cg05951221, cg21566642
plot(METH['cg01940273',], METH['cg05951221',], xlim = c(0,1), ylim = c(0,1), pch = 19, col = alpha('black', 0.5))
plot(METH['cg01940273',], METH['cg21566642',], xlim = c(0,1), ylim = c(0,1), pch = 19, col = alpha('black', 0.5))
plot(METH['cg05951221',], METH['cg21566642',], xlim = c(0,1), ylim = c(0,1), pch = 19, col = alpha('black', 0.5))

setwd("<not_shown>")
tiff(filename = 'correlogram.tiff', res = 300, width = 4, height = 4.1, units = 'in')
corrplot(cor(t(METH), method = 'spearman'), method="circle", order = 'hclust', type = 'lower', diag = F, 
         tl.cex = 0.5, cl.cex = 0.5)
dev.off()

####################### Sex #######################

X = as.data.frame(t(METH))
X$sex = phenotype[match(colnames(METH), phenotype$ERF_ID),]$ERF_gender
X$smoke = phenotype[match(colnames(METH), phenotype$ERF_ID),]$smoking_3cat
df = melt(X, id.vars = c('sex', 'smoke'))
df$smoke = factor(df$smoke, c('non_smoker', 'past_smoker', 'current_smoker'))

p = ggplot(data = na.omit(df), aes(x = sex, y = value, fill = smoke)) + geom_boxplot(alpha = 0.8) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_wrap(~variable) +
  scale_fill_manual(values = c('bisque1', 'gray82', 'dimgray')) + xlab("")

setwd("<not_shown>")
ggsave(filename = 'sex.tiff', plot = p, device = 'tiff', width = 6.7, height = 6.2, units = 'in', dpi = 300)

pvals = sapply(1:nrow(METH), function(x) summary(lm(X[,x] ~ sex + smoke, data = X))$coefficients['sexM','Pr(>|t|)'])
names(pvals) = rownames(METH)
p.adjust(pvals, 'bonferroni')
# cg01940273  cg03636183  cg05575921  cg05951221  cg06126421  cg09935388  cg12803068  cg12876356  cg13039251  cg15693572 
# 1.000000000 1.000000000 1.000000000 1.000000000 0.137838223 0.003732179 1.000000000 1.000000000 1.000000000 1.000000000 
# cg21566642  cg22132788  cg23576855 
# 0.239765620 0.402675209 1.000000000 

####################### Pack years, time since cessation #######################

X = as.data.frame(t(METH))
X$smoke = phenotype[match(colnames(METH), phenotype$ERF_ID),]$smoking_3cat
X$n_cigarettes = phenotype[match(colnames(METH), phenotype$ERF_ID),]$n_cigarettes_current
X$time_since_cessation = phenotype[match(colnames(METH), phenotype$ERF_ID),]$age - phenotype[match(colnames(METH), phenotype$ERF_ID),]$age_quit_cigarettes
X$n_cigarettes[X$smoke != 'current_smoker'] = 0
X$time_since_cessation[X$smoke != 'past_smoker'] = NA
df = melt(X, id.vars = c('n_cigarettes', 'smoke', 'time_since_cessation'))
df$smoke = factor(df$smoke, c('non_smoker', 'past_smoker', 'current_smoker'))

table(is.na(df$n_cigarettes), df$smoke)
#        non_smoker past_smoker current_smoker
# FALSE        923         949           1170
# TRUE           0           0             52

p = ggplot(data = df[!is.na(df$n_cigarettes) & !is.na(df$smoke),], aes(x = n_cigarettes, y = value, col = smoke)) + geom_point(alpha = 0.8) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_wrap(~variable, ncol = 6) +
  scale_color_manual(values = c('bisque1', 'gray82', 'dimgray')) + xlab("") + geom_smooth(method = 'lm', linetype = 2, col = 'red2')
setwd("<not_shown>")
ggsave(filename = 'ncigaretteaday.tiff', plot = p, device = 'tiff', width = 7, height = 4, units = 'in', dpi = 300)

pvals = sapply(1:nrow(METH), function(x) summary(lm(X[!is.na(X$n_cigarettes) & !is.na(df$smoke),x] ~ X$n_cigarettes[!is.na(X$n_cigarettes) & !is.na(df$smoke)]))$coefficients[2,'Pr(>|t|)'])
names(pvals) = colnames(X)[1:13]

p = ggplot(data = df[!is.na(df$time_since_cessation) & !is.na(df$smoke),], aes(x = time_since_cessation, y = value, col = smoke)) + geom_point(alpha = 0.8) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_wrap(~variable, ncol = 6) +
  scale_color_manual(values = c('gray82', 'dimgray')) + xlab("") + geom_smooth(method = 'lm', linetype = 2, col = 'red2')
setwd("<not_shown>")
ggsave(filename = 'timesincecessation.tiff', plot = p, device = 'tiff', width = 7, height = 4, units = 'in', dpi = 300)

pvals = sapply(1:nrow(METH), function(x) summary(lm(X[!is.na(X$time_since_cessation) & !is.na(df$smoke),x] ~ X$time_since_cessation[!is.na(X$time_since_cessation) & !is.na(df$smoke)]))$coefficients[2,'Pr(>|t|)'])
names(pvals) = colnames(X)[1:13]

length(X[!is.na(X$n_cigarettes) & !is.na(X$smoke), 1]) # 237
length(X[!is.na(X$time_since_cessation) & !is.na(X$smoke), 1]) # 72

####################### Prediction modelling #######################

# Remove technical replicates
rep = names(which(table(colnames(METH_corrected)) > 1))
indices = sapply(1:length(rep), function(x) which(colnames(METH) == rep[x]))
REP = METH[, c(indices[1,], indices[2,])]
REP_corrected = METH_corrected[, c(indices[1,], indices[2,])]
X = METH[,-indices[2,]]
X_corrected = METH_corrected[,-indices[2,]]

# Remove NAs in phenotype
smoke3cat = phenotype[match(colnames(X), phenotype$ERF_ID),]$smoking_3cat
X = X[, smoke3cat != '']
X_corrected = X_corrected[, smoke3cat != '']
smoke3cat = smoke3cat[smoke3cat != '']
  
# Extract phenotypes
smoke2cat = smoke3cat
smoke2cat[smoke2cat == ''] = NA
smoke2cat[smoke2cat %in% c('non_smoker', 'past_smoker')] = 'past_or_non_smoker'
smoke2cat = factor(smoke2cat, c('past_or_non_smoker', 'current_smoker'))
table(smoke2cat, exclude = F)
# smoke2cat
# past_or_non_smoker     current_smoker
#   142                 90             

smoke3cat = factor(smoke3cat, c('non_smoker', 'past_smoker', 'current_smoker'))
table(smoke3cat, exclude = F)
# smoke3cat
# non_smoker    past_smoker current_smoker
#   71             71             90      

setwd("<not_shown>")
UMtools::export_bigmat(bigmat = cbind(t(round(X, 5)), phenotype[match(colnames(X), phenotype$ERF_ID),]), 'dataset_uncorrected.txt')

setwd("<not_shown>")
UMtools::export_bigmat(bigmat = cbind(t(round(X_corrected, 5)), phenotype[match(colnames(X), phenotype$ERF_ID),]), 'dataset_corrected.txt')

################################ Predict uncorrected ################################

############## 2-cat
X = data.frame(t(X))
pred2cat = predict(M2cat, X, type = "response")
pred2cat = as.numeric(pred2cat > 0.5)
table(pred2cat, smoke2cat, exclude = F)
# pred2cat past_or_non_smoker current_smoker
# 0                126             36
# 1                 16             54

confusionMatrix(factor(pred2cat), factor(as.numeric(smoke2cat == 'current_smoker')), mode = "everything", positive="1")
# F1 : 0.6750

############## 3-cat
pred3cat = predict(M3cat, X, type = "probs") 
pred3cat = unlist(sapply(1:nrow(pred3cat), function(x) unname(which.max(pred3cat[x,]))))

table(pred3cat, smoke3cat, exclude = F)
#          smoke3cat
# pred3cat non_smoker past_smoker current_smoker
# 1         58          26             15       
# 2          6          35             19       
# 3          7          10             56       

tmp = smoke3cat; levels(tmp) = c(1, 2, 3)
F1_Score_macro_weighted(tmp, pred3cat) # 0.6487559

################################ Predict corrected ################################

############## 2-cat
X = data.frame(t(X_corrected))
pred2cat = predict(M2cat, X, type = "response")
pred2cat = as.numeric(pred2cat > 0.5)
table(pred2cat, smoke2cat, exclude = F)
# pred2cat past_or_non_smoker current_smoker
# 0                139             66
# 1                  3             24

confusionMatrix(factor(pred2cat), factor(as.numeric(smoke2cat == 'current_smoker')), mode = "everything", positive="1")
# F1 : 0.4103

############## 3-cat
pred3cat = predict(M3cat, X, type = "probs") 
pred3cat = unlist(sapply(1:nrow(pred3cat), function(x) unname(which.max(pred3cat[x,]))))

table(pred3cat, smoke3cat, exclude = F)
# pred3cat non_smoker past_smoker current_smoker
# 1         65          66             52
# 2          5           3              9
# 3          1           2             29   

tmp = smoke3cat; levels(tmp) = c(1, 2, 3)
F1_Score_macro_weighted(tmp, pred3cat) # 0.4612271

################################ Microarray Vs seq ################################

# Extract data
markers = c('cg12876356', 'cg21566642', 'cg09935388', 'cg15693572', 'cg05951221', 'cg03636183',
            'cg01940273', 'cg13039251', 'cg23576855', 'cg06126421', 'cg12803068', 'cg22132788', 'cg05575921')
df = M2cat$data[,c(markers, 'Smoking')]
mod = glm(formula = Smoking ~ ., data = df, family = 'binomial')
plot(coef(mod), coef(M2cat)); abline(0,1, lty = 2) # Same model

table(as.numeric(predict(M2cat, newX = df, type = 'response') > 0.5), df$Smoking)
#      0    1
# 0 3174  212
# 1   79  299

confusionMatrix(data = table(as.numeric(predict(M2cat, newX = df, type = 'response') > 0.5), df$Smoking), mode = 'everything')$byClass['F1']
#        F1 
# 0.9561681 

# Effect sizes - Odds ratio from individual models
effects = sapply(1:13, function(x) standardize_parameters(glm(formula = Smoking ~ ., data = df[,c(x, 14)], family = 'binomial'), exp = T)$Std_Odds_Ratio[2])
names(effects) = markers
sort(effects)
# cg05575921 cg21566642 cg01940273 cg05951221 cg03636183 cg23576855 cg06126421 cg09935388 cg12876356 cg15693572 cg12803068 
# 0.2176432  0.2573151  0.3661238  0.3805879  0.3982737  0.5157057  0.5328573  0.5550271  0.7057990  1.6928272  1.8635150 
# cg13039251 cg22132788 
# 1.9273534  2.0470204 

df$Smoking = factor(df$Smoking)
levels(df$Smoking) = c('non-smoker', 'past or current smoker')
df = melt(df, id.vars = "Smoking")
df$variable = as.character(df$variable)
df$variable = factor(df$variable, levels = LEV)
p = ggplot(data = na.omit(df), aes(x = variable, y = value, fill = Smoking)) + geom_boxplot() + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = c('bisque1', 'dimgray')) + xlab("")
setwd("<not_shown>")
ggsave(filename = '2cat_silvana_450K.tiff', plot = p, device = 'tiff', width = 6.7, height = 6.2, units = 'in', dpi = 300)

# UNCORRECTED
df_2 = as.data.frame(t(METH))
y = phenotype[match(colnames(METH), phenotype$ERF_ID),]$smoking_3cat
y[y == ''] = NA; y[y %in% c('current_smoker', 'past_smoker')] = 'past or current smoker'; y[y == 'non_smoker'] = 'non-smoker'
df_2$smoke = factor(y)
df_2 = melt(df_2, id.vars = "smoke")
df_2$variable = as.character(df_2$variable)
df_2$variable = factor(df_2$variable, levels = LEV)
df$platform = 'seq'
df_2$platform = 'microarray'
names(df_2)[1] = 'Smoking'
df_all = rbind(df, df_2)
df_all$label = paste(df_all$platform, df_all$Smoking, sep = ': ')
df_all$label = factor(df_all$label, levels = c('microarray: non-smoker', 'seq: non-smoker',
                                               'microarray: past or current smoker', 'seq: past or current smoker'))
p = ggplot(data = na.omit(df_all), aes(x = value, y = label, fill = label)) + geom_density_ridges() + facet_wrap(~variable) + 
  theme_minimal() + scale_fill_manual(values = c('slategray1', 'slategray4', 'papayawhip', 'lightsteelblue')) + xlab('') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position= 'none')  + xlim(0,1)
setwd("<not_shown>")
ggsave(filename = '450K_sequncorr_comp.tiff', plot = p, device = 'tiff', width = 6.7, height = 6.2, units = 'in', dpi = 300)

# CORRECTED
df_2 = as.data.frame(t(METH_corrected))
y = phenotype[match(colnames(METH_corrected), phenotype$ERF_ID),]$smoking_3cat
y[y == ''] = NA; y[y %in% c('current_smoker', 'past_smoker')] = 'past or current smoker'; y[y == 'non_smoker'] = 'non-smoker'
df_2$smoke = factor(y)
df_2 = melt(df_2, id.vars = "smoke")
df_2$variable = as.character(df_2$variable)
df_2$variable = factor(df_2$variable, levels = LEV)
df$platform = 'seq'
df_2$platform = 'microarray'
colnames(df_2)[1] = 'Smoking'
df_all = rbind(df, df_2)
df_all$label = paste(df_all$platform, df_all$Smoking, sep = ': ')
df_all$label = factor(df_all$label, levels = c('microarray: non-smoker', 'seq: non-smoker',
                                               'microarray: past or current smoker', 'seq: past or current smoker'))
p = ggplot(data = na.omit(df_all), aes(x = value, y = label, fill = label)) + geom_density_ridges() + facet_wrap(~variable) + 
  theme_minimal() + scale_fill_manual(values = c('slategray1', 'slategray4', 'papayawhip', 'lightsteelblue')) + xlab('') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position= 'none') + xlim(0,1)
setwd("<not_shown>")
ggsave(filename = '450K_seqcorr_comp.tiff', plot = p, device = 'tiff', width = 6.7, height = 6.2, units = 'in', dpi = 300)

################################ Joint Microarray-seq model (2-cat) - uncorrected ################################

# Build joint dataset
markers = c('cg12876356', 'cg21566642', 'cg09935388', 'cg15693572', 'cg05951221', 'cg03636183',
            'cg01940273', 'cg13039251', 'cg23576855', 'cg06126421', 'cg12803068', 'cg22132788', 'cg05575921')
df = M2cat$data[,c(markers, 'Smoking')]
df$technology = 'microarray'
df$Smoking = as.factor(df$Smoking)
levels(df$Smoking) = c('past_or_non_smoker', 'current_smoker')
df2 = as.data.frame(t(X)[,markers])
df2$Smoking = smoke2cat
df2$technology = 'sequencing'
df_global = rbind(df, df2)
df_global$technology = as.factor(df_global$technology)

## mQTL
smoke = df$Smoking
smoke = factor(as.numeric(smoke == 'current_smoker'))
levels(smoke) = c('bisque1', 'dimgray'); smoke = as.character(smoke)
plot(df[,'cg05575921'], df[,'cg23576855'], xlim = c(0,1), ylim = c(0,1), col = alpha(smoke, 0.3), pch = 19)
##

set.seed(123)
seeds <- vector(mode = "list", length = 11)
for(i in 1:10) seeds[[i]] <- sample.int(1000, 1)
seeds[[11]] <- sample.int(1000, 1)
train_control <- trainControl(method="repeatedcv", number=5, repeats = 2, savePredictions = TRUE, classProbs = TRUE, summaryFunction = f1, seeds = seeds)
model <- train(Smoking ~ technology + technology:cg12876356 + technology:cg21566642 + technology:cg09935388 + technology:cg15693572 +
               technology:cg05951221 + technology:cg03636183 + technology:cg01940273 + technology:cg13039251 + technology:cg23576855 +
               technology:cg06126421 + technology:cg12803068 + technology:cg22132788 + technology:cg05575921,
               data=df_global, trControl=train_control, method="glm", family = "binomial", metric = 'F1')
PREDICTIONS = model$pred
PREDICTIONS$technology = df_global$technology[PREDICTIONS$rowIndex]
CONF_MAT = table(PREDICTIONS$pred, PREDICTIONS$obs, PREDICTIONS$Resample, PREDICTIONS$technology)
CONF_MAT = lapply(1:length(unique(PREDICTIONS$Resample)), function(x) lapply(1:2, function(y) CONF_MAT[,,x,y]))
PERFORMANCE = sapply(1:length(unique(PREDICTIONS$Resample)), function(x) sapply(1:2, function(y) confusionMatrix(data = CONF_MAT[[x]][[y]], mode = 'everything')$byClass['F1']))
rownames(PERFORMANCE) = c('microarray', 'seq')
colnames(PERFORMANCE) = as.character(sapply(1:5, function(x) paste(x, 1:2, sep = '_')))
#                  1_1       1_2       2_1       2_2       3_1       3_2       4_1       4_2       5_1       5_2
# microarray 0.9567854 0.9581465 0.9509434 0.9574944 0.9581465 0.9548585 0.9601803 0.9554044 0.9540317 0.9527382
# seq        0.8055556 0.7719298 0.7936508 0.8000000 0.8064516 0.8354430 0.8235294 0.8307692 0.8214286 0.8000000

round(rowMeans(PERFORMANCE), 3)
# microarray   seq 
# 0.956      0.809 
round(rowSds(PERFORMANCE), 3)
# [1] 0.003 0.019

# Final model
summary(model$finalModel)
# Coefficients:
# Estimate Std. Error z value Pr(>|z|)
# (Intercept)                      -2.89464    1.50117  -1.928 0.053823 .
# technologysequencing              5.05522    6.01492   0.840 0.400658
# technologymicroarray:cg12876356   4.57264    1.54346   2.963 0.003051 **
# technologysequencing:cg12876356  -2.55343    7.30752  -0.349 0.726771
# technologymicroarray:cg21566642  -4.27201    1.57303  -2.716 0.006612 **
# technologysequencing:cg21566642  -2.33539    3.15073  -0.741 0.458559
# technologymicroarray:cg09935388  -3.03172    1.39250  -2.177 0.029468 *
# technologysequencing:cg09935388  -0.15680    2.76528  -0.057 0.954781
# technologymicroarray:cg15693572   2.12076    0.70061   3.027 0.002470 **
# technologysequencing:cg15693572   0.46380    1.65669   0.280 0.779513
# technologymicroarray:cg05951221   4.89178    1.81891   2.689 0.007158 **
# technologysequencing:cg05951221  -0.04798    3.25069  -0.015 0.988223
# technologymicroarray:cg03636183   6.43007    1.55836   4.126 3.69e-05 ***
# technologysequencing:cg03636183  -0.45527    3.24688  -0.140 0.888488
# technologymicroarray:cg01940273  -5.46031    2.14948  -2.540 0.011076 *
# technologysequencing:cg01940273   0.84024    4.64855   0.181 0.856562
# technologymicroarray:cg13039251   5.61982    1.28518   4.373 1.23e-05 ***
# technologysequencing:cg13039251   0.08810    2.14564   0.041 0.967248
# technologymicroarray:cg23576855  -1.63228    0.48362  -3.375 0.000738 ***
# technologysequencing:cg23576855  -0.55051    1.34731  -0.409 0.682833
# technologymicroarray:cg06126421   4.78682    1.09257   4.381 1.18e-05 ***
# technologysequencing:cg06126421   5.39726    2.65677   2.032 0.042203 *
# technologymicroarray:cg12803068  -9.37791    1.51869  -6.175 6.62e-10 ***
# technologysequencing:cg12803068  -0.41507    3.31247  -0.125 0.900282
# technologymicroarray:cg22132788  13.44069    2.14658   6.261 3.81e-10 ***
# technologysequencing:cg22132788   2.51454    3.07876   0.817 0.414077
# technologymicroarray:cg05575921 -18.52626    1.17214 -15.805  < 2e-16 ***
# technologysequencing:cg05575921  -6.75439    2.59116  -2.607 0.009142 **

setwd("<not_shown>")
M2cat_MPS = model$finalModel
save(M2cat_MPS, file='Vidaki_et_al_2cat.rda')

# ROC curves - 2 categories
y = as.numeric(df_global$Smoking == 'current_smoker')[df_global$technology == 'sequencing']
pred <- prediction(predictions = predict(model$finalModel, type = "response")[df_global$technology == 'sequencing'], 
                   labels = y)
roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
acc.perf = performance(pred, measure = "acc"); plot(acc.perf)
auc_ROCR <- performance(pred, measure = "auc"); auc_ROCR@y.values[[1]] # 0.8362285
plot(roc.perf, main = paste("AUC =", round(auc_ROCR@y.values[[1]], 5))); abline(a = 0, b = 1, lty = 2)

y_pred = as.numeric(predict(model$finalModel, type = 'response') > 0.5)
y_real = as.numeric(df_global$Smoking == 'current_smoker')
setwd("<not_shown>")
tiff(filename = 'microarray_2cat_conf.tiff', res = 300, width = 4, height = 4.1, units = 'in')
gplots::balloonplot(table(y_real[df_global$technology == 'microarray'], 
                          y_pred[df_global$technology == 'microarray']), main = '')
dev.off()

tiff(filename = 'seq_2cat_conf.tiff', res = 300, width = 4, height = 4.1, units = 'in')
gplots::balloonplot(table(y_real[df_global$technology == 'sequencing'], 
                          y_pred[df_global$technology == 'sequencing']), main = '')
dev.off()

df_rep = data.frame(t(REP))
df_rep$technology = 'microarray'
Prob = data.frame(predict(model, df_rep, type = "prob"))
cat = factor(sapply(1:nrow(Prob), function(x) which.max(Prob[x,])), 1:2)
tiff(filename = 'reps_2cat_conf.tiff', res = 300, width = 4, height = 4.1, units = 'in')
gplots::balloonplot(table(unname(cat[1:6]), unname(cat[7:12])), main = '')
dev.off()

###
y_pred__2cat = y_pred[df_global$technology == 'sequencing']
y_real__2cat = y_real[df_global$technology == 'sequencing']
##

################################ Joint Microarray-seq model (3-cat) - uncorrected ################################

markers = c('cg12876356', 'cg21566642', 'cg09935388', 'cg15693572', 'cg05951221', 'cg03636183',
            'cg01940273', 'cg13039251', 'cg23576855', 'cg06126421', 'cg12803068', 'cg22132788', 'cg05575921')
df = M2cat$data[,c(markers, 'Smoking')]
u = log(as.numeric(M3cat$fitted.values[,2]/M3cat$fitted.values[,3]))
w = as.numeric(cbind(1, as.matrix(df[,colnames(coef(M3cat))[-1]])) %*% (coef(M3cat)[1,]-coef(M3cat)[2,]))
positions = sapply(1:length(u), function(x) which.min((u[x]-w)^2))
df = df[positions,]
df$Smoking = sapply(1:nrow(M3cat$fitted.values), function(x) which((M3cat$fitted.values+M3cat$residuals)[x,] == 1))
df$Smoking = factor(df$Smoking)
table(df$Smoking)
# 1    2     3 
# 1243 1332  364
table(df$Smoking, predict(M3cat))
#     1   2   3
# 1 969 272   2
# 2 407 868  57
# 3  24  97 243
F1_Score_macro_weighted(df$Smoking, predict(M3cat)) # 0.7090622

mod = multinom(formula = Smoking ~ ., data = df)
plot(coef(mod), coef(M3cat)[,colnames(coef(mod))]); abline(0, 1, lty = 2)
df$technology = 'microarray'
levels(df$Smoking) = c('non_smoker', 'past_smoker', 'current_smoker')
df2 = as.data.frame(t(X)[,markers])
df2$Smoking = smoke3cat
df2$technology = 'sequencing'
df_global = rbind(df, df2)
df_global$technology = as.factor(df_global$technology)

set.seed(124)
seeds <- vector(mode = "list", length = 11)
for(i in 1:10) seeds[[i]] <- sample.int(1000, 3)
seeds[[11]] <- sample.int(1000, 1)
train_control <- trainControl(method="repeatedcv", number=5, repeats = 2, savePredictions = TRUE, classProbs = TRUE, summaryFunction = macrof1, seeds = seeds)
model <- train(Smoking ~ technology + technology:cg12876356 + technology:cg21566642 + technology:cg09935388 + technology:cg15693572 +
                 technology:cg05951221 + technology:cg03636183 + technology:cg01940273 + technology:cg13039251 + technology:cg23576855 +
                 technology:cg06126421 + technology:cg12803068 + technology:cg22132788 + technology:cg05575921,
               data=df_global, trControl=train_control, method="multinom", metric = 'F1')
PREDICTIONS = model$pred
PREDICTIONS$technology = df_global$technology[PREDICTIONS$rowIndex]
PERFORMANCE = rbind(sapply(unique(PREDICTIONS$Resample), function(x) macrof1(PREDICTIONS[PREDICTIONS$Resample == x & 
                                                                             PREDICTIONS$technology == 'microarray',])),
                    sapply(unique(PREDICTIONS$Resample), function(x) macrof1(PREDICTIONS[PREDICTIONS$Resample == x & 
                                                                             PREDICTIONS$technology == 'sequencing',])))
rownames(PERFORMANCE) = c('microarray', 'seq')
colnames(PERFORMANCE) = as.character(sapply(1:5, function(x) paste(x, 1:2, sep = '_')))
#                  1_1       1_2       2_1       2_2       3_1       3_2       4_1       4_2       5_1       5_2
# microarray 0.7195563 0.6984020 0.6876356 0.6922730 0.7216709 0.6877306 0.7360022 0.7021016 0.6602734 0.7061514
# seq        0.6680097 0.5564081 0.6056434 0.6132393 0.6513932 0.6450975 0.6580164 0.6910317 0.5274814 0.4829053

round(rowMeans(PERFORMANCE), 3)
# microarray        seq 
# 0.701      0.610 
round(rowSds(PERFORMANCE), 3)
# [1] 0.021 0.068

# Final model
model$finalModel
# Coefficients:
#                 (Intercept) technologysequencing `technologymicroarray:cg12876356` `technologysequencing:cg12876356`
# past_smoker       17.02616           -13.722483                          3.412244                          1.317331
# current_smoker    11.70163            -6.655559                          8.128004                         -1.711276
#                 `technologymicroarray:cg21566642` `technologysequencing:cg21566642` `technologymicroarray:cg09935388`
# past_smoker                            -11.81367                         -2.200274                        -0.7533385
# current_smoker                         -13.32560                         -3.872069                        -3.8956487
#                 `technologysequencing:cg09935388` `technologymicroarray:cg15693572` `technologysequencing:cg15693572`
# past_smoker                            -1.542352                        -0.5654644                         -3.459018
# current_smoker                         -1.213950                         2.6430503                         -1.671068
#                 `technologymicroarray:cg05951221` `technologysequencing:cg05951221` `technologymicroarray:cg03636183`
# past_smoker                            0.1173349                         -7.827782                         -4.286820
# current_smoker                         5.0676710                         -3.898821                          3.311741
#                 `technologysequencing:cg03636183` `technologymicroarray:cg01940273` `technologysequencing:cg01940273`
# past_smoker                           -0.4725785                          4.387662                          3.811822
# current_smoker                        -0.9408307                         -3.623739                          3.432423
#                 `technologymicroarray:cg13039251` `technologysequencing:cg13039251` `technologymicroarray:cg23576855`
# past_smoker                           -0.6488779                         0.6326889                        -0.8771121
# current_smoker                         3.3439474                         0.5651880                        -3.0449158
#                 `technologysequencing:cg23576855` `technologymicroarray:cg06126421` `technologysequencing:cg06126421`
# past_smoker                            2.0516838                        -7.5499663                         -3.784824
# current_smoker                         0.3010023                        -0.1801977                          3.083727
#                 `technologymicroarray:cg12803068` `technologysequencing:cg12803068` `technologymicroarray:cg22132788`
# past_smoker                             6.880009                          6.005643                         -5.886301
# current_smoker                         -4.179963                          2.839890                         12.035594
#               `technologysequencing:cg22132788` `technologymicroarray:cg05575921` `technologysequencing:cg05575921`
# past_smoker                           -5.5264838                         -7.276351                         0.3950404
# current_smoker                        -0.5994861                        -27.724727                        -6.1742781

setwd("<not_shown>")
M3cat_MPS = model$finalModel
save(M3cat_MPS, file='Vidaki_et_al_3cat.rda')

y_pred = predict(model$finalModel)
levels(y_pred) = c(1,2,3)
y_real = df_global$Smoking
levels(y_real) = c(1,2,3)

###
y_pred__3cat = y_pred[df_global$technology == 'sequencing']
y_real__3cat = y_real[df_global$technology == 'sequencing']
##

setwd("<not_shown>")
tiff(filename = 'microarray_3cat_conf.tiff', res = 300, width = 4, height = 4.1, units = 'in')
gplots::balloonplot(table(y_real[df_global$technology == 'microarray'], 
                          y_pred[df_global$technology == 'microarray']), main = '')
dev.off()

tiff(filename = 'seq_3cat_conf.tiff', res = 300, width = 4, height = 4.1, units = 'in')
gplots::balloonplot(table(y_real[df_global$technology == 'sequencing'], 
                          y_pred[df_global$technology == 'sequencing']), main = '')
dev.off()
Prob = data.frame(predict(model, df_global[df_global$technology == 'microarray',], type = "prob"))
Prob$reported = df_global$Smoking[df_global$technology == 'microarray']
p = ggtern(data=Prob, aes(x=non_smoker,y=past_smoker, z=current_smoker, col = reported)) +  geom_point(alpha = 0.1, size = 0.5)+
  limit_tern(1.1,1.1,1.1) + scale_color_manual(values = c('forestgreen', 'goldenrod1', 'darkslateblue'))
setwd("<not_shown>")
ggsave(filename = 'tern_microarray.tiff', plot = p, 
       device = 'tiff', width = 6.7/1.5, height = 6/2, units = 'in', dpi = 300)

Prob = data.frame(predict(model, df_global[df_global$technology == 'sequencing',], type = "prob"))
Prob$reported = df_global$Smoking[df_global$technology == 'sequencing']
p = ggtern(data=Prob, aes(x=non_smoker,y=past_smoker, z=current_smoker, col = reported)) +  geom_point(alpha = 0.5, size = 1)+
  limit_tern(1.1,1.1,1.1) + scale_color_manual(values = c('forestgreen', 'goldenrod1', 'darkslateblue'))

setwd("<not_shown>")
ggsave(filename = 'tern_seq.tiff', plot = p, 
       device = 'tiff', width = 6.7/1.5, height = 6/2, units = 'in', dpi = 300)

# REP
df_rep = data.frame(t(REP))
df_rep$technology = 'microarray'
Prob = data.frame(predict(model, df_rep, type = "prob"))
cat = factor(sapply(1:nrow(Prob), function(x) which.max(Prob[x,])), 1:3)
tiff(filename = 'reps_3cat_conf.tiff', res = 300, width = 4, height = 4.1, units = 'in')
gplots::balloonplot(table(unname(cat[1:6]), unname(cat[7:12])), main = '')
dev.off()

table(y_real__2cat, y_real__3cat)
#               y_real__3cat
# y_real__2cat  1  2  3
# 0            71 71  0
# 1             0  0 90

matrices = lapply(1:3, function(x) table(y_pred__2cat, y_pred__3cat, y_real__3cat)[,,x])
setwd("<not_shown>")
tiff(filename = '2cat_vs_3cat_1', res = 300, width = 4, height = 4.1, units = 'in')
gplots::balloonplot(matrices[[1]], main = '')
dev.off()
tiff(filename = '2cat_vs_3cat_2', res = 300, width = 4, height = 4.1, units = 'in')
gplots::balloonplot(matrices[[2]], main = '')
dev.off()
tiff(filename = '2cat_vs_3cat_3', res = 300, width = 4, height = 4.1, units = 'in')
gplots::balloonplot(matrices[[3]], main = '')
dev.off()

################################ Joint Microarray-seq model (2-cat) - corrected ################################

# Build joint dataset
markers = c('cg12876356', 'cg21566642', 'cg09935388', 'cg15693572', 'cg05951221', 'cg03636183',
            'cg01940273', 'cg13039251', 'cg23576855', 'cg06126421', 'cg12803068', 'cg22132788', 'cg05575921')
df = M2cat$data[,c(markers, 'Smoking')]
df$technology = 'microarray'
df$Smoking = as.factor(df$Smoking)
levels(df$Smoking) = c('past_or_non_smoker', 'current_smoker')
df2 = as.data.frame(t(X_corrected)[,markers])
df2$Smoking = smoke2cat
df2$technology = 'sequencing'
df_global = rbind(df, df2)
df_global$technology = as.factor(df_global$technology)

set.seed(125)
seeds <- vector(mode = "list", length = 11)
for(i in 1:10) seeds[[i]] <- sample.int(1000, 1)
seeds[[11]] <- sample.int(1000, 1)
train_control <- trainControl(method="repeatedcv", number=5, repeats = 2, savePredictions = TRUE, classProbs = TRUE, summaryFunction = f1, seeds = seeds)
model <- train(Smoking ~ technology + technology:cg12876356 + technology:cg21566642 + technology:cg09935388 + technology:cg15693572 +
                 technology:cg05951221 + technology:cg03636183 + technology:cg01940273 + technology:cg13039251 + technology:cg23576855 +
                 technology:cg06126421 + technology:cg12803068 + technology:cg22132788 + technology:cg05575921,
               data=df_global, trControl=train_control, method="glm", family = "binomial", metric = 'F1')

PREDICTIONS = model$pred
PREDICTIONS$technology = df_global$technology[PREDICTIONS$rowIndex]
CONF_MAT = table(PREDICTIONS$pred, PREDICTIONS$obs, PREDICTIONS$Resample, PREDICTIONS$technology)
CONF_MAT = lapply(1:length(unique(PREDICTIONS$Resample)), function(x) lapply(1:2, function(y) CONF_MAT[,,x,y]))
PERFORMANCE = sapply(1:length(unique(PREDICTIONS$Resample)), function(x) sapply(1:2, function(y) confusionMatrix(data = CONF_MAT[[x]][[y]], mode = 'everything')$byClass['F1']))
rownames(PERFORMANCE) = c('microarray', 'seq')
colnames(PERFORMANCE) = as.character(sapply(1:5, function(x) paste(x, 1:2, sep = '_')))
#                  1_1       1_2       2_1       2_2       3_1       3_2       4_1       4_2       5_1       5_2
# microarray 0.9517345 0.9493294 0.9561338 0.9539623 0.9630189 0.9571751 0.9530303 0.9631856 0.9538229 0.9571865
# seq        0.7272727 0.7547170 0.8750000 0.8181818 0.8571429 0.7500000 0.8307692 0.8571429 0.7812500 0.7792208

round(rowMeans(PERFORMANCE), 3)
# microarray        seq 
# 0.956      0.803 
round(rowSds(PERFORMANCE), 3)
# [1] 0.005 0.052

################################ Joint Microarray-seq model (3-cat) - corrected ################################

markers = c('cg12876356', 'cg21566642', 'cg09935388', 'cg15693572', 'cg05951221', 'cg03636183',
            'cg01940273', 'cg13039251', 'cg23576855', 'cg06126421', 'cg12803068', 'cg22132788', 'cg05575921')
df = M2cat$data[,c(markers, 'Smoking')]
u = log(as.numeric(M3cat$fitted.values[,2]/M3cat$fitted.values[,3]))
w = as.numeric(cbind(1, as.matrix(df[,colnames(coef(M3cat))[-1]])) %*% (coef(M3cat)[1,]-coef(M3cat)[2,]))
positions = sapply(1:length(u), function(x) which.min((u[x]-w)^2))
df = df[positions,]
df$Smoking = sapply(1:nrow(M3cat$fitted.values), function(x) which((M3cat$fitted.values+M3cat$residuals)[x,] == 1))
df$Smoking = factor(df$Smoking)
df$technology = 'microarray'
levels(df$Smoking) = c('non_smoker', 'past_smoker', 'current_smoker')
df2 = as.data.frame(t(X_corrected)[,markers])
df2$Smoking = smoke3cat
df2$technology = 'sequencing'
df_global = rbind(df, df2)
df_global$technology = as.factor(df_global$technology)

set.seed(126)
seeds <- vector(mode = "list", length = 11)
for(i in 1:10) seeds[[i]] <- sample.int(1000, 3)
seeds[[11]] <- sample.int(1000, 1)
train_control <- trainControl(method="repeatedcv", number=5, repeats = 2, savePredictions = TRUE, classProbs = TRUE, summaryFunction = macrof1, seeds = seeds)
model <- train(Smoking ~ technology + technology:cg12876356 + technology:cg21566642 + technology:cg09935388 + technology:cg15693572 +
                 technology:cg05951221 + technology:cg03636183 + technology:cg01940273 + technology:cg13039251 + technology:cg23576855 +
                 technology:cg06126421 + technology:cg12803068 + technology:cg22132788 + technology:cg05575921,
               data=df_global, trControl=train_control, method="multinom", metric = 'F1')

PREDICTIONS = model$pred
PREDICTIONS$technology = df_global$technology[PREDICTIONS$rowIndex]
PERFORMANCE = rbind(sapply(unique(PREDICTIONS$Resample), function(x) macrof1(PREDICTIONS[PREDICTIONS$Resample == x & 
                                                                                           PREDICTIONS$technology == 'microarray',])),
                    sapply(unique(PREDICTIONS$Resample), function(x) macrof1(PREDICTIONS[PREDICTIONS$Resample == x & 
                                                                                           PREDICTIONS$technology == 'sequencing',])))
rownames(PERFORMANCE) = c('microarray', 'seq')
colnames(PERFORMANCE) = as.character(sapply(1:5, function(x) paste(x, 1:2, sep = '_')))
#                  1_1       1_2       2_1       2_2       3_1       3_2       4_1       4_2       5_1       5_2
# microarray 0.7047757 0.6932766 0.7112812 0.6936548 0.7007619 0.7147112 0.7061721 0.7038863 0.7065953 0.6834256
# seq        0.7500300 0.5367425 0.5579194 0.4881635 0.6408350 0.6464702 0.6992941 0.5914653 0.6050328 0.5734495

round(rowMeans(PERFORMANCE), 3)
# microarray        seq 
# 0.702      0.609 
round(rowSds(PERFORMANCE), 3)
# [1] 0.009 0.078

################################ Joint Microarray-seq model Ordinal logistic regression (3-cat) - uncorrected ################################

markers = c('cg12876356', 'cg21566642', 'cg09935388', 'cg15693572', 'cg05951221', 'cg03636183',
            'cg01940273', 'cg13039251', 'cg23576855', 'cg06126421', 'cg12803068', 'cg22132788', 'cg05575921')
df = M2cat$data[,c(markers, 'Smoking')]
u = log(as.numeric(M3cat$fitted.values[,2]/M3cat$fitted.values[,3]))
w = as.numeric(cbind(1, as.matrix(df[,colnames(coef(M3cat))[-1]])) %*% (coef(M3cat)[1,]-coef(M3cat)[2,]))
positions = sapply(1:length(u), function(x) which.min((u[x]-w)^2))
df = df[positions,]
df$Smoking = sapply(1:nrow(M3cat$fitted.values), function(x) which((M3cat$fitted.values+M3cat$residuals)[x,] == 1))
df$Smoking = factor(df$Smoking)
table(df$Smoking)
# 1    2     3 
# 1243 1332  364
table(df$Smoking, predict(M3cat))
#     1   2   3
# 1 969 272   2
# 2 407 868  57
# 3  24  97 243
df$technology = 'microarray'
levels(df$Smoking) = c('non_smoker', 'past_smoker', 'current_smoker')
df2 = as.data.frame(t(X)[,markers])
df2$Smoking = smoke3cat
df2$technology = 'sequencing'
df_global = rbind(df, df2)
df_global$technology = as.factor(df_global$technology)
df_global$Smoking = factor(df_global$Smoking, ordered = T)

set.seed(127)
seeds <- vector(mode = "list", length = 11)
for(i in 1:10) seeds[[i]] <- sample.int(1000, 5)
seeds[[11]] <- sample.int(1000, 1)
train_control <- trainControl(method="repeatedcv", number=5, repeats = 2, savePredictions = TRUE, classProbs = TRUE, summaryFunction = macrof1, seeds = seeds)
model <- train(Smoking ~ technology + technology:cg12876356 + technology:cg21566642 + technology:cg09935388 + technology:cg15693572 +
                 technology:cg05951221 + technology:cg03636183 + technology:cg01940273 + technology:cg13039251 + technology:cg23576855 +
                 technology:cg06126421 + technology:cg12803068 + technology:cg22132788 + technology:cg05575921,
               data=df_global, trControl=train_control, method="polr", metric = 'F1')
PREDICTIONS = model$pred
PREDICTIONS$technology = df_global$technology[PREDICTIONS$rowIndex]

PERFORMANCE = rbind(sapply(unique(PREDICTIONS$Resample), function(x) macrof1(PREDICTIONS[PREDICTIONS$Resample == x & 
                                                                                           PREDICTIONS$technology == 'microarray',])),
                    sapply(unique(PREDICTIONS$Resample), function(x) macrof1(PREDICTIONS[PREDICTIONS$Resample == x & 
                                                                                           PREDICTIONS$technology == 'sequencing',])))
rownames(PERFORMANCE) = c('microarray', 'seq')
colnames(PERFORMANCE) = as.character(sapply(1:5, function(x) paste(x, 1:2, sep = '_')))
#                  1_1       1_2       2_1       2_2       3_1       3_2       4_1       4_2       5_1       5_2
# microarray 0.6730685 0.6542443 0.6530468 0.6592651 0.7010594 0.6936557 0.6424644 0.7042681 0.6479169 0.6404519
# seq        0.5597321 0.5373220 0.5273975 0.5709948 0.6460865 0.6398450 0.5217272 0.5150730 0.5733528 0.5656705

round(rowMeans(PERFORMANCE), 3)
# microarray        seq 
# 0.667      0.566 
round(rowSds(PERFORMANCE), 3)
# [1] 0.024 0.046

# Final model
model$finalModel

############################################################################################

sessionInfo()
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.6 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=nl_NL.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=nl_NL.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=nl_NL.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=nl_NL.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1 UMtools_0.1                                        dbscan_1.1-10                                     
# [4] Sushi_1.34.0                                       biomaRt_2.52.0                                     zoo_1.8-10                                        
# [7] EMCluster_0.2-14                                   Matrix_1.4-1                                       minfi_1.42.0                                      
# [10] bumphunter_1.38.0                                  locfit_1.5-9.6                                     iterators_1.0.14                                  
# [13] foreach_1.5.2                                      Biostrings_2.64.0                                  XVector_0.36.0                                    
# [16] SummarizedExperiment_1.26.1                        Biobase_2.56.0                                     MatrixGenerics_1.8.1                              
# [19] matrixStats_0.62.0                                 GenomicRanges_1.48.0                               GenomeInfoDb_1.32.3                               
# [22] IRanges_2.30.0                                     S4Vectors_0.34.0                                   BiocGenerics_0.42.0                               
# [25] effectsize_0.8.1                                   MLmetrics_1.1.3                                    HandTill2001_1.0.1                                
# [28] ROCR_1.0-11                                        minpack.lm_1.2-2                                   caret_6.0-93                                      
# [31] lattice_0.20-45                                    nnet_7.3-17                                        MASS_7.3-58.1                                     
# [34] corrplot_0.92                                      scales_1.2.0                                       RColorBrewer_1.1-3                                
# [37] ggtern_3.3.5                                       ggridges_0.5.4                                     ggplot2_3.3.6                                     
# [40] gplots_3.1.3                                       data.table_1.14.2                                 
# 
# loaded via a namespace (and not attached):
# [1] utf8_1.2.2                proto_1.0.0               tidyselect_1.2.0          RSQLite_2.2.15            AnnotationDbi_1.58.0      grid_4.2.1               
# [7] BiocParallel_1.30.3       pROC_1.18.0               munsell_0.5.0             codetools_0.2-18          preprocessCore_1.58.0     future_1.28.0            
# [13] withr_2.5.0               colorspace_2.0-3          filelock_1.0.2            rstudioapi_0.13           robustbase_0.95-0         bayesm_3.1-4             
# [19] listenv_0.8.0             GenomeInfoDbData_1.2.8    datawizard_0.6.3          bit64_4.0.5               rhdf5_2.40.0              parallelly_1.32.1        
# [25] vctrs_0.4.1               generics_0.1.3            ipred_0.9-13              BiocFileCache_2.4.0       R6_2.5.1                  illuminaio_0.38.0        
# [31] bitops_1.0-7              rhdf5filters_1.8.0        cachem_1.0.6              reshape_0.8.9             DelayedArray_0.22.0       assertthat_0.2.1         
# [37] BiocIO_1.6.0              gtable_0.3.0              globals_0.16.1            timeDate_4021.106         rlang_1.0.6               genefilter_1.78.0        
# [43] splines_4.2.1             rtracklayer_1.56.1        ModelMetrics_1.2.2.2      GEOquery_2.64.2           yaml_2.3.5                reshape2_1.4.4           
# [49] GenomicFeatures_1.48.3    tensorA_0.36.2            tools_4.2.1               lava_1.7.0                nor1mix_1.3-0             ellipsis_0.3.2           
# [55] siggenes_1.70.0           latex2exp_0.9.5           Rcpp_1.0.9                plyr_1.8.7                sparseMatrixStats_1.8.0   progress_1.2.2           
# [61] zlibbioc_1.42.0           purrr_0.3.4               RCurl_1.98-1.8            prettyunits_1.1.1         rpart_4.1.16              openssl_2.0.2            
# [67] magrittr_2.0.3            hms_1.1.1                 xtable_1.8-4              XML_3.99-0.10             mclust_5.4.10             gridExtra_2.3            
# [73] compiler_4.2.1            tibble_3.1.8              KernSmooth_2.23-20        crayon_1.5.1              tzdb_0.3.0                tidyr_1.2.0              
# [79] lubridate_1.8.0           DBI_1.1.3                 dbplyr_2.2.1              rappdirs_0.3.3            compositions_2.0-4        readr_2.1.2              
# [85] cli_3.4.1                 quadprog_1.5-8            insight_0.18.6            gower_1.0.0               pkgconfig_2.0.3           GenomicAlignments_1.32.1 
# [91] recipes_1.0.2             xml2_1.3.3                annotate_1.74.0           hardhat_1.2.0             rngtools_1.5.2            multtest_2.52.0          
# [97] beanplot_1.3.1            prodlim_2019.11.13        doRNG_1.8.2               scrime_1.3.5              stringr_1.4.0             digest_0.6.29            
# [103] parameters_0.19.0         base64_2.0                DelayedMatrixStats_1.18.0 restfulr_0.0.15           curl_4.3.2                Rsamtools_2.12.0         
# [109] gtools_3.9.3              rjson_0.2.21              lifecycle_1.0.3           nlme_3.1-159              Rhdf5lib_1.18.2           askpass_1.1              
# [115] limma_3.52.2              fansi_1.0.3               pillar_1.8.0              KEGGREST_1.36.3           fastmap_1.1.0             httr_1.4.3               
# [121] DEoptimR_1.0-11           survival_3.4-0            glue_1.6.2                bayestestR_0.13.0         png_0.1-7                 bit_4.0.4                
# [127] class_7.3-20              stringi_1.7.8             HDF5Array_1.24.2          blob_1.2.3                caTools_1.18.2            memoise_2.0.1            
# [133] dplyr_1.0.9               future.apply_1.9.1      
