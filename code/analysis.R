

library(ggplot2)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(GGally)


# Arjune S. Dhanekula, MD in Dave Marcinek's lab had some RNA seq done on
# mice aorta (destructive sampling) 
# old and young
# treated with ELAM or control

# data from Azenta Life Sci company were aquired on Illumina technology.
# Azenta had done data analysis 

# I will take the data normalized by Azenta (assuming they remove techncial variation that I might be unaware of or unable to handle) and work forward with the goals of 
# 1. characterization of the data/samples (Azenta found differential genes among pairwise samples (ie, old v. young))
# 2. compare to senescence-associated gene expression profiles
# 3. comapre to any available clock models of the aorta
# 4. other stuff?

rm(list=ls())

getwd() 
setwd() 

# Arjune sent meta data for the samples:
samps <- read.csv("RNAseq_Mice_info.csv", head=T) # I simplified the meta data into variables
head(samps)

# Arjune found a mouse with a misspecified sex [his email]: "-One of my mice was mistakenly labeled as a M when it was a F. It’s only one mouse, but I thought you may want to update your data to still ensure that there are no sex differences (GP115-5, Age ELAM, should be a F and was mistakenly labeled as M)"

samps[samps$sample=='GP115-5', ]
samps$sex[samps$sample=='GP115-5'] <- 'F'


dir()
dir('AzenaFiles/DEG/deseq2')

# code to join data:
azentaComparisons <- print(dir('AzenaFiles/DEG/deseq2'))
azentaComparisons <- azentaComparisons[azentaComparisons != "DESeq2_column_definitions.pdf"]
dataPathI <- 'AzenaFiles/DEG/deseq2/' 
dataPathII <- '/counts/rlog_transformed_counts.csv'

dList <- list()
for(i in 1:length(azentaComparisons)) {
dList[[i]] <- read.csv(file=paste0(dataPathI, azentaComparisons[i], dataPathII) , head=T, row.names=1) }
names(dList) <- LETTERS[1:5]

head(dList[[1]])
genes <- sapply(dList, function(x) rownames(x))
genes <- Reduce(union, genes) 

d <- cbind(dList[[1]][genes, !colnames(dList[[1]])%in% colnames(dList[[2]])], dList[[2]][genes, !colnames(dList[[2]])%in% colnames(dList[[1]])])
head(d)
colSums(is.na(d))

d <- cbind(d, dList[[3]][genes, !colnames(dList[[3]]) %in% colnames(d)])
colSums(is.na(d))

d <- cbind(d, dList[[4]][genes, !colnames(dList[[4]]) %in% colnames(d)])
colSums(is.na(d))

d <- cbind(d, dList[[5]][genes, !colnames(dList[[5]]) %in% colnames(d)])
colSums(is.na(d))

head(d)
d <- as.data.frame(t(d))
d[1:4,1:4]

s <- samps$sample
s[1:10]
rownames(d)

s2 <- gsub('-', '.', s)
s2[!s2 %in% rownames(d)] # these samples aquired an additional tag from Aventa, likely needed to be re-run, modify the names to match, and merge the data:

rownames(d) <- gsub('resub.', '', rownames(d))
samps$sample <- gsub('-', '.', s)

head(samps)
samps$sample %in% rownames(d)
d$sample <- rownames(d)

d <- merge(samps, d)

table(genes %in% colnames(d))
genes <- genes[genes %in% colnames(d)]

d$ageGroup <- factor(ifelse(d$age >10, 'old', 'young'), levels=c('young', 'old'))

d[1:4,1:10]
d <- d[ ,colnames(d)!='cage.1']

table(grepl('NA', colnames(d)))

d[ ,grepl('NA', colnames(d))]
table(colSums(is.na(d)), 'NA in colname'=grepl('NA', colnames(d)))
d <- d[ ,!grepl('NA', colnames(d))]

d$cage <- as.factor(d$cage)

save(d, genes, file='merged.data')


ggpairs(d[ ,c('age', 'ageGroup', 'sex', 'treat', 'cage')], aes(color=treat))



rm(list=ls())
################################################
load('merged.data')

hist(colSums(is.na(d[ ,genes])), xlab='samples with missing values', ylab='genes', main='') 

completeGenes <- genes[colSums(is.na(d[ ,genes]))==0] # 20,607 geens with no missing data, of 21,326 genes

par(mfrow=c(1,2))
boxplot(d[ ,completeGenes[1:100]], xlab='genes', ylab='normalized reads from first 1000 genes', cex=0.5, main='genes are normalized')
boxplot(t(d[ ,completeGenes]), xlab='samples', ylab='normalized reads from 20,607 complete genes', cex=0.5, main='sample-wide bias apparently removed')


table(d$age, d$treat)

pca <- prcomp(d[ ,completeGenes ], scale=T)

tw <- AssocTests::tw(pca$sdev, eigenL = nrow(d))
sigEigs <- tw$SigntEigenL # 3 signficant eigenvalues

scree <- data.frame('prop variance' = pca$sdev^2/sum(pca$sdev^2), PC=as.factor(1:nrow(d)))
head(scree)

summary(pca) # the fisrt 3 PCs explain 58.7% of the variance

sum(scree$prop.variance[1:3])

ggplot(scree, aes(y=prop.variance, x=PC, fill=PC))+
  geom_bar(stat='identity')+
  theme_bw(base_size = 14)+
  scale_fill_manual(values = c("1" = "#ff8d33", '2'='#ff8d33', '3'='#ff8d33'))+
  theme(, legend.position='none')+
  ylab('variance explained (prop)')

pca <- pca$x
head(pca)
PCs <- paste0('PC', 1:6)

head(d[ ,1:10])
pca <- data.frame(treat=as.factor(d$treat), sex=as.factor(d$sex), age=d$age, cage=as.factor(d$cage), pca)        
head(pca)

pca$ageGroup <- factor(ifelse(pca$age >10, 'old', 'young'), levels=c('young', 'old'))

ggpairs(pca[ ,c('age', 'ageGroup', 'sex', 'treat', 'cage', PCs)], aes(color=treat)) # PC1 ~ age group x treatment, nice!

table(d$cage) # consider using a random effects model on cage (although many cages are n=1, there may be enough of the cages with n=2, 3 to fit intercepts to each)

library(lme4)
library(lmerTest)

PCs <- paste0('PC', 1:sigEigs)
p <- list()

for(i in 1:3) {
mm <- lmer(pca[ ,PCs[i]] ~ sex + treat * ageGroup + (1|cage), data=pca)
a <- car::Anova(mm)
p[[i]] <- a$`Pr(>Chisq)` }
names(p) <- PCs

p <- do.call(rbind, p)
colnames(p) <- rownames(a)

p <- apply(p, 2, function(x) p.adjust(x, 'bonferroni'))

p 
# PC1 and PC3 sig for treat, 
# PCs 1 and 2 for age group
# PC1, PC3* (* barely) for age group x treatment


pca$condition <- paste( pca$ageGroup, pca$treat)

pca$condition <- factor(pca$condition, levels=c("young control", "young ELAM", "old control", "old ELAM"))

pca$condition

PCs <- paste0('PC', 1:3)
l <- pivot_longer(pca, cols=PCs, values_to = 'eigenvalue', names_to = 'PC')
head(l)

freeScalePCs <- ggplot(l, aes(eigenvalue, x=condition, color=condition, fill=condition))+
  geom_point()+
  geom_boxplot(alpha=0.5)+
  theme_bw(base_size = 16)+
  scale_fill_manual(values=c('#ff005e', '#93003a', '#93c4d2', '#00429d'))+
  scale_color_manual(values=c('#ff005e', '#93003a', '#93c4d2', '#00429d'))+
  facet_wrap(~PC, nrow=1, scales='free')+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())

fixedScalePCs <- ggplot(l, aes(eigenvalue, x=condition, color=condition, fill=condition))+
  geom_point()+
  geom_boxplot(alpha=0.5)+
  theme_bw(base_size = 16)+
  scale_fill_manual(values=c('#ff005e', '#93003a', '#93c4d2', '#00429d'))+
  scale_color_manual(values=c('#ff005e', '#93003a', '#93c4d2', '#00429d'))+
  facet_wrap(~PC, nrow=1)+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())

freeScalePCs
fixedScalePCs
ggarrange(freeScalePCs, fixedScalePCs, ncol=1, common.legend = T, legend='right')


######################################
## pairwise correlation among samples
######################################
cMat <- cor(t(d[ ,completeGenes]))
rownames(cMat) <- as.character(1:24)
colnames(cMat) <- as.character(1:24)
pheatmap(cMat)

ann <- data.frame(ageGroup=d$ageGroup, treatment=d$treat, sex=d$sex)
head(ann)
str(ann)
ann_colors = list(ageGroup = c(`young`='lightgreen', `old`='darkgreen'), treatment=c(ELAM='purple', control='lightblue'), sex=c(F='orange', M='red'))
ann
rownames(ann) <- colnames(cMat) 

pheatmap(cMat, scale='none', annotation_col = ann, annotation_colors = ann_colors, main='', fontsize=6, fontsize_row=1, fontsize_col = 1)
# similarity among samples looks inflated, in part due to co-linarity among mzs.  The PCA is likely better-reflecting the relationship among the samples in the metabolome




###########################################################################
# univariate analysis:
###########################################################################

# consider a mixed model, to account for the 1 to 3 measurements of mice that share a cage
mean(table(d$cage))
range(table(d$cage))
table(table(d$cage)) # the distribution of sample sizes among the cages. 8 single mice per cage, 5 cages with 2 mice, 2 cages with three mice

## Here is a useful discussion of the issue/choice to assign a random effect to a variable with n=1 to 3 samples:
# from: https://stats.stackexchange.com/questions/242821/how-will-random-effects-with-only-1-observation-affect-a-generalized-linear-mixe
#In general, you have an issue with identifiability. Linear models with a random effect assigned to a parameter with only one measurement can't distinguish between the random effect and the residual error. A typical linear mixed effect equation will look like:
  
#  𝐸=𝛽+𝜂𝑖+𝜖𝑗

# Where 𝛽is the fixed effect, 𝜂𝑖 is the random effect for level 𝑖, and 𝜖𝑗is the residual variability for the 𝑗th measurement. When you have only one observation of a level with a random effect, it is difficult to distinguish between 𝜂and 𝜖. You will (typically) be fitting a variance or standard deviation to 𝜂 and 𝜖, so with only one measure per individual, you will be not be as certain that you have an accurate estimate for 𝑆𝐷(𝜂) and 𝑆𝐷(𝜖), but the estimate of the sum of the variances (𝑣𝑎𝑟(𝜂)+𝑣𝑎𝑟(𝜖)) should be relatively robust.

# On to the practical answer: If you have about 1/3 of your observations with a single observation per individual, you are probably OK overall. The rest of the population should provide a reasonable estimate for 𝑆𝐷(𝜂) and 𝑆𝐷(𝜖) , and these individuals should be minor contributors overall. On the other hand, if you have all individuals at a specific fixed effect and random effect with a single measure (e.g. for your example, perhaps all of a population-- perhaps that means species for you), then you would trust the result less.

## NOTE: tried models like that below with a random effect of cage, but, not surprisingly, it failed to converge, likely due to all of the single instances of many of the cages. (n = 1 to 3 sampled mice per cage).  As the abive comment points out, failure to converge is not unexpected, and may be ok.  I compared the result (using i=1) with AIC:

aicList <- list()

cores <- detectCores() 
cl <- makeCluster(cores)
registerDoParallel(cl)

for(i in 1:40) {
m <- lm(d[ ,completeGenes[i]] ~ sex + treat * ageGroup, data=d) 
mm <- lmer(d[ ,completeGenes[i]] ~ sex + treat * ageGroup + (1|cage), data=d) 
aicList[[i]] <- AIC(m, mm) } # a much lower AIC for the linear model without the random effect. 

aic <- t(do.call(cbind, aicList))
aic <- aic[rownames(aic)=='AIC', ]
plot(aic, main='AIC')
abline(0,1) # the aic is systematically lower for the linear model w/o the random effect of cage
# this isn't a surprise considering the discussion above - takeaway, adding a cage effect does not improve the information attributed to the model parameters, esp. compared to the cost of the additional parameter
######################################################################################################


######################################################################################################
# linear univariate model:
betas <-list()
p <- list()

cores <- detectCores() 
cl <- makeCluster(cores)
registerDoParallel(cl)

for(i in 1:length(completeGenes)) {
  m <- lm(d[ ,completeGenes[i]] ~ sex + treat * ageGroup, data=d) 
  s <- summary(m)
  betas[[i]] <- s$coefficients[-c(1), 1]
  a <- car::Anova(m)
  p[[i]] <- a$`Pr(>F)` }

names(p) <- completeGenes
names(betas) <- completeGenes
betaNames <- rownames(s$coefficients[-c(1), ])


save(p, a, d, completeGenes, betas, betaNames, file='univariateResults')
#######################################################




rm(list=ls())
#######################################################
load('univariateResults')

p <- do.call(rbind, p)
colnames(p) <- rownames(a)

head(p)
fdr <- as.data.frame(apply(p, 2, function(x) p.adjust(x, 'fdr')))
colSums(fdr<=0.05) # number of genes with effects of each model term

x <- completeGenes[fdr$`treat:ageGroup`<=0.05]
x <- t(d[ ,x])

dev.off()
pheatmap(x, scale='row', fontsize = 1)
colnames(x) <- rownames(d)

ann <- data.frame(ageGroup=d$ageGroup, treatment=d$treat, sex=d$sex)
head(ann)
str(ann)
ann_colors = list(ageGroup = c(`young`='lightgreen', `old`='darkgreen'), treatment=c(ELAM='purple', control='lightblue'), sex=c(M='pink', `F`='red'))
ann
rownames(ann) <- colnames(x) 
pheatmap(x, scale='row', annotation_col = ann, annotation_colors = ann_colors, main='', fontsize=6, fontsize_row=1, fontsize_col = 1)


colSums(fdr<=0.05) # number of genes with effects of each model term

table(fdr$treat<=0.05, fdr$`treat:ageGroup`<=0.05)
fisher.test(table(fdr$treat<=0.05, fdr$`treat:ageGroup`<=0.05)) # the genes with treatment effects are enriched among those with treatment by age effects.
plot(table(fdr$treat<=0.05, fdr$`treat:ageGroup`<=0.05))



table(table(d$cage))



rm(list=ls())
################################################################################################
### look at the relationship between the genes responsive to ELAM, or ELAM by age, and the senescence-associated SenMayo gene set
################################################################################################

################################################
load('merged.data')
hist(colSums(is.na(d[ ,genes])), xlab='samples with missing values', ylab='genes', main='') 

completeGenes <- genes[colSums(is.na(d[ ,genes]))==0] # 20,607 geens with no missing data, of 21,326 genes
d <- d[ ,!colnames(d) %in% genes[!genes%in%completeGenes]] # remove incomplete genes from data

# downloaded the SenMayo gene set at:
# https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/SAUL_SEN_MAYO.html
s <- read.csv('sen_mayo_genes.csv', head=T) # 125 (mouse?) genes, which this site says map to 117 mouse genes.
head(s)
s$Gene.Id <- as.character(s$Gene.Id)

completeGenes[1]

library(stringr)
converts <- data.frame('gene'=completeGenes, 'Gene.Id'=as.character(as.numeric(str_sub(completeGenes, start= -6))))

head(converts)
s$Gene.Id %in% converts$Gene.Id

write.table(completeGenes, file='ENSMUSGs_to_convert.csv', quote=F, row.names = F, sep=",",  col.names=FALSE) # use external tool (https://biit.cs.ut.ee/gprofiler/convert) to convert ENSMUSG to ENSP, which also grabs gene symbols

# load gene ID and gene name info:
ids <- read.csv('gProfiler_mmusculus_11-19-2024_9-04-22 AM.csv', head=T) 
head(ids)

s$Symbol %in% ids$name
s[!s$Symbol %in% ids$name, ]

table('Number_Match' =s$Gene.Id %in% converts$Gene.Id, 'Symbol_Match' =s$Symbol %in% ids$name) # of 125 sen mayo genes, 22 do not match either gene symbols or gene numbers.  NOTE: there are 8 genes from SenMayo that didn't map to a mouse gene to begin with (https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/SAUL_SEN_MAYO.html). So this is more like of 117 genes, 14 do not map

s[!s$Gene.Id %in% converts$Gene.Id & !s$Symbol %in% ids$name, ] # these are the genes that do not match
allArjuneGenes <- as.character(as.numeric(str_sub(genes, start= -6)))

table('Number_Match' =s$Gene.Id %in% converts$Gene.Id, 'Symbol_Match' =s$Symbol %in% ids$name) 
table('Number_Match' =s$Gene.Id %in% allArjuneGenes, 'Symbol_Match' =s$Symbol %in% ids$name) 

table(s$Symbol=='', is.na(s$Gene.Id)) # 7 genes in sen mayo set are NA


## take what I can get of the sen mayo genes and analyze?  # I may miss genes in the Arjune set that are in sen mayo under a different name...
table('Number_Match' =s$Gene.Id %in% converts$Gene.Id, 'Symbol_Match' =s$Symbol %in% ids$name)# 86+13 = 99 matches by gene symbol

head(ids)
head(converts)
rownames(ids) <- ids$initial_alias
rownames(converts) <- converts$gene

ids$sen_mayo <- as.factor(ifelse(converts$Gene.Id %in% s$Gene.Id | ids$name %in% s$Symbol, 'sen', 'bkd'))
table(ids$sen_mayo) # 114 sen mayo genes of 117 (97%), with 20493 background genes

dir()
load('univariateResults')

p <- do.call(rbind, p)
colnames(p) <- rownames(a)

head(p)
fdr <- as.data.frame(apply(p, 2, function(z) p.adjust(z, 'fdr')))
colSums(fdr<=0.05)

head(ids)
head(fdr)
fdr$initial_alias <- rownames(fdr)

o <- merge(fdr, ids)
head(o)

colnames(o)[2:5] <- paste0(colnames(o)[2:5], '_FDR')
save(o, file='genesForGO')


##########################################################################
# plot to depict the effect size and direction (betas) among genes, and in relation to SenMayo genes
##########################################################################
betas[[1]]
b <- as.data.frame(do.call(rbind, betas))
head(b)
colnames(b) <- paste0(colnames(b), '_beta')

head(o)
rownames(o) <- o$initial_alias
head(b)
b <- b[rownames(o), ]
o <- as.data.frame(cbind(o, b))


head(o)

o <- o[ ,colnames(o) != 'namespace']
o <- o[ ,colnames(o) != 'Residuals']
write.table(o, quote=F, row.names = F, file='table.for.Arjune.txt', sep='\t')

table(o$`treatELAM:ageGroupold_beta`[o$`treat:ageGroup_FDR`<=0.05]>0, o$sen_mayo[o$`treat:ageGroup_FDR`<=0.05])


levels(o$sen_mayo) <- c('background', 'senescence gene')

o$hit <- as.factor(o$`treat:ageGroup_FDR`<=0.5)
levels(o$hit) <- c('NS', 'FDR<5%')

ggplot(o, aes(x=-`treatELAM:ageGroupold_beta`, y=-log(`treat:ageGroup_FDR`,10), color=sen_mayo, label=name))+
  geom_point(aes(shape=hit))+
  theme_bw(base_size = 18)+
  scale_shape_manual(values=c(1,19))+
  facet_wrap(~sen_mayo)+
  ylab('-log10 FDR')+
  xlab('effect of ELAM on old-age difference')+
  labs(shape='age x ELAM effect') +
  labs(color='gene set') +
  ggrepel::geom_text_repel(max.overlaps = 8)


table('age-assoc'=o$ageGroup_FDR<=0.05, o$sen_mayo, 'up with age'=o$ageGroupold_beta<0)
## of the 56-age associated SenMayo genes, the expression of 53 increase with age

table(o$hit,  o$sen_mayo, 'age x treat'=o$`treat:ageGroup_FDR`<=0.05)

table(o$sen_mayo, o$`treat:ageGroup_FDR`<=0.05, 'up'=o$`treatELAM:ageGroupold_beta`<0) 
## of the 43-age associated SenMayo genes with treatment-depdent effects, all of them increase with age !!!


dev.off() # clear the large plot if you're done with it.
##########################################################################

table(o$sen_mayo)



##########################################################################
# intersection among SneMayo genes and genes with model effects:
##########################################################################
table(o$sen_mayo, o$`treat:ageGroup_FDR`<=0.05) # 43 genes in SenMayo are among the treatment x age effect genes
table(o$sen_mayo, o$treat_FDR<=0.05) # 31 among the treatment genes
table(o$sen_mayo, o$ageGroup_FDR<=0.05) # 56 among the age genes

fisher.test(table(o$sen_mayo, o$`treat:ageGroup_FDR`<=0.05)) # sen mayo genes are enriched among the treatment-dependent age-effect genes
fisher.test(table(o$sen_mayo, o$treat_FDR<=0.05))
fisher.test(table(o$sen_mayo, o$ageGroup_FDR <=0.05)) # sen mayo genes are enriched amogn the age-associated genes

table(o$sen_mayo)/nrow(o) # change of being a sen mayo gene = 0.005532
table(o$ageGroup_FDR<=0.05)/nrow(o) # chance of being an age-assoc gene = 0.295
both <- 0.294 * 0.005532101 # change of being both
both*nrow(o) # 33 expected

table(o$`treat:ageGroup_FDR`<=0.05)/nrow(o) # change of being a treatment x age gene = 0.186
both <- 0.2105595 * 0.005532101 # change of being both
both*nrow(o) # 24 expected number of genes that would be both treat x age and sen mayo

table(o$treat_FDR<=0.05)/nrow(o) # chance of being a treatment gene = 0.305
both <- 0.305 * 0.005532101 # change of being both
both*nrow(o) # 35 expected number of genes that would be both treat x age and sen mayo


plot(table(o$sen_mayo, 'ELAM x age (FDR<5%)'=o$`treat:ageGroup_FDR`<=0.05), main='genes intersection between senMayo and ELAM x age', cex=1.5) # visual representation of the significant over-represetnation of senMayo genes among the genes assocaited with ELAM x age

##########################################################################



##########################################################################
## heatmap
##########################################################################
senMayo <- o$sen_mayo # make a named vector of sen mayo genes
names(senMayo) <- o$initial_alias
table(senMayo)


x <- completeGenes[fdr$`treat:ageGroup`<=0.05]
x <- t(d[ ,x])
pheatmap(x, scale='row', fontsize = 1)
colnames(x) <- rownames(d)

ann <- data.frame(ageGroup=d$ageGroup, treatment=d$treat, sex=d$sex)
head(ann)
str(ann)

row.ann <- data.frame(sen_mayo = senMayo[rownames(x)])
head(row.ann)

row.ann

ann_colors = list(ageGroup = c(`young`='lightgreen', `old`='darkgreen'), treatment=c(ELAM='purple', control='lightblue'), sen_mayo=c(`senescence gene`='black', background='white'))
 

pheatmap(x, scale='row', annotation_col = ann, annotation_colors = ann_colors, main='', fontsize=12, fontsize_row=1, fontsize_col = 1, annotation_row = row.ann, clustering_method= 'complete')


##########################################################################


##########################################################################
## long form data frame for transcriptome-wide modeling
l <- pivot_longer(d, cols=all_of(completeGenes), values_to='expression', names_to='gene')
head(l)

l$sen_mayo <- senMayo[l$gene]
table(l$sen_mayo)
levels(l$sen_mayo) <- c('bkd', 'sen')

m <- aggregate(expression ~ gene + ageGroup + treat, l, mean)
head(m)

m$sen_mayo <- senMayo[m$gene]
head(m)

m$sen_mayo

ggplot(m, aes(y=expression, x=ageGroup, group=gene))+
  geom_line(alpha=0.1)+
  theme_classic()+
  ggh4x::facet_grid2(~treat ~sen_mayo, scales = "free_y", independent = "y") # wow!


head(m)

m$sen_mayo <- factor(m$sen_mayo, levels=c('senescence gene', 'background'))
levels(m$sen_mayo) <- c('SenMayo', 'background')
m$treat[m$treat=='control'] <- 'NT'
m$treat <- factor(m$treat, levels=c('NT', 'ELAM'))


ggplot(m, aes(y=expression, x=treat, fill=sen_mayo))+
  geom_line(alpha=0.2, aes(group=gene))+
  geom_violin(alpha=0.5, linewidth = 0, scale='area')+
  theme_classic(base_size = 16)+
  labs(fill='gene set')+
  ggh4x::facet_grid2(~sen_mayo ~ageGroup , scales = "free_y", independent = "y")



ggplot(m, aes(y=expression, x=sen_mayo, fill=treat))+
  geom_line(alpha=0.2, aes(group=gene))+
  geom_violin(alpha=0.5, linewidth = 0, scale='area')+
  theme_classic(base_size = 16)+
  labs(fill='gene set')+
  ggh4x::facet_grid2(~treat ~ageGroup , scales = "free_y", independent = "y")


bkdA <- ggplot(subset(m, sen_mayo=='background'), aes(y=expression, x=treat, fill=sen_mayo))+
  geom_line(alpha=0.2, aes(group=gene))+
  geom_violin(alpha=0.5, linewidth = 0, scale='area')+
  theme_classic(base_size = 16)+
  theme(legend.position = 'none')+
  facet_wrap(~ageGroup)

bkdB <- ggplot(subset(m, sen_mayo=='background'), aes(y=expression, x=ageGroup, fill=sen_mayo))+
  geom_line(alpha=0.2, aes(group=gene))+
  geom_violin(alpha=0.5, linewidth = 0, scale='area')+
  theme_classic(base_size = 16)+
  theme(legend.position = 'none')+
  facet_wrap(~treat)

ggarrange(bkdA, bkdB, common.legend = T, legend='none')




head(m)

w <- pivot_wider(m, names_from = treat, values_from = expression)
head(w)

w$DiffELAM <- w$ELAM-w$control
head(w)

ggplot(w, aes(DiffELAM, fill = sen_mayo, colour = sen_mayo)) + 
  facet_wrap(~ageGroup, scales='free', ncol=1)+
  geom_histogram(aes(y=..density..), alpha=0.7, position="identity")+
  scale_color_manual(values=c('SenMayo'='red', 'background'='grey'))+
  scale_fill_manual(values=c('SenMayo'='red', 'background'='grey'))+
  theme_classic(base_size = 16)+
  xlim(c(-2, 1.2))

levels(w$ageGroup)
levels(w$sen_mayo)
w$sen_mayo <-  relevel(w$sen_mayo, ref='background')
summary(lm(DiffELAM ~ sen_mayo * ageGroup, w))





ggplot(subset(m, sen_mayo=='senescence gene'), aes(y=expression, x=treat, fill=sen_mayo))+
  geom_violin(alpha=0.2, linewidth = 0)+
  geom_line(alpha=0.5, aes(group=gene))+
  theme_classic(base_size = 16)+
  theme(legend.position = 'none') +
  ylab('sen mayo gene expression')+
  xlab('')+
  scale_fill_manual(values=c('purple'))+
  facet_wrap(~ageGroup)







head(l)


s <- summary(lmer(expression ~ ageGroup * treat * sen_mayo + sex + (1|gene) + (1|cage), data=l))
par(mar=c(10,20,10,10))
barplot(s$coefficients[-c(1),1], horiz=T, las=1)
par(mar=c(5,4,4,2))
s

summary(lmer(expression ~ ageGroup * sen_mayo + sex + (1|gene) + (1|cage), data=l[l$treat=='ELAM', ]))

p2

a <- car::Anova(lmer(expression ~ ageGroup * treat * sen_mayo + sex + (1|gene) + (1|cage), data=l))
a

p1 <- ggplot(m, aes(y=expression, x=sen_mayo))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=0.5, width=0.1)+
  theme_classic()

ggpubr::ggarrange(p1, p2, p3, nrow=1, widths=c(1, 2, 2))

levels(m$sen_mayo) <- c('background gene', 'senescence gene')

p3 <- ggplot(m, aes(y=expression, x=ageGroup, group=gene, color=sen_mayo))+
  geom_line(alpha=0.5)+
  theme_classic()+
  ggh4x::facet_grid2(~sen_mayo ~treat , scales = "free_y", independent = "y")
p3 + theme_classic(base_size = 16) + xlab('')


senPlot <- ggplot(subset(m, sen_mayo=='senescence gene'), aes(y=expression, x=ageGroup, group=gene, color=sen_mayo))+
  geom_line(alpha=0.5)+
  theme_classic(base_size = 16)+
  theme(legend.position = 'none') +
  xlab('')+
  ylab('sen mayo gene expression')+
  scale_color_manual(values=c('red'))+
  facet_wrap(~treat)

bkdPlot <- ggplot(subset(m, sen_mayo=='background gene'), aes(y=expression, x=ageGroup, group=gene, color=sen_mayo))+
  geom_line(alpha=0.1)+
  theme_classic(base_size = 16)+
  theme(legend.position = 'none') +
  ylab('background gene expression')+
  xlab('')+
  scale_color_manual(values=c('black'))+
  facet_wrap(~treat)


ggarrange(senPlot, bkdPlot, ncol=1)


# plot to show difference btw control and treatment at each age:
senTA <- ggplot(subset(m, sen_mayo=='senescence gene'), aes(y=expression, x=treat, group=gene, color=sen_mayo))+
  geom_line(alpha=0.5)+
  theme_classic(base_size = 16)+
  theme(legend.position = 'none') +
  ylab('sen mayo gene expression')+
  xlab('')+
  scale_color_manual(values=c('purple'))+
  facet_wrap(~ageGroup)

bkdTA <- ggplot(subset(m, sen_mayo=='background gene'), aes(y=expression, x=treat, group=gene, color=sen_mayo))+
  geom_line(alpha=0.1)+
  theme_classic(base_size = 16)+
  theme(legend.position = 'none') +
  ylab('background gene expression')+
  xlab('')+
  scale_color_manual(values=c('black'))+
  facet_wrap(~ageGroup)


ggarrange(senTA, bkdTA, ncol=1)



senPlot <- ggplot(subset(m, sen_mayo=='senescence gene'), aes(y=expression, x=ageGroup, fill=sen_mayo))+
  geom_violin(alpha=0.2, linewidth = 0)+
  geom_line(alpha=0.5, aes(group=gene))+
  theme_classic(base_size = 16)+
  theme(legend.position = 'none') +
  xlab('')+
  ylab('sen mayo gene expression')+
  scale_fill_manual(values=c('red'))+
  facet_wrap(~treat)

senTA <- ggplot(subset(m, sen_mayo=='senescence gene'), aes(y=expression, x=treat, fill=sen_mayo))+
  geom_violin(alpha=0.2, linewidth = 0)+
  geom_line(alpha=0.5, aes(group=gene))+
  theme_classic(base_size = 16)+
  theme(legend.position = 'none') +
  ylab('sen mayo gene expression')+
  xlab('')+
  scale_fill_manual(values=c('purple'))+
  facet_wrap(~ageGroup)


ggarrange(senPlot, senTA, nrow=1)
