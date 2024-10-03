library(ashr)
library(mashr)
library("dplyr")
library("tidyr")
library(data.table)
library("reshape2")
library("matrixStats")


setwd("/storage/home/snagpal3/Decanalization/Amplification/")

# field
fd="GeneticSex"
# select exposure pair
sel="GeneticSex-Testosterone"


fielddict <- read.csv(paste0("/storage/home/snagpal3/Decanalization/Amplification/Fields_categories_subcategories_withGroupCode_newCode_color.csv"), header =T)
freq <- fread(paste0("/storage/home/snagpal3/Decanalization/Data/UKB_freq_WB.txt"), header = T)

# example of mash for exposure genetic sex and testosterone
print(sel)
female_df <- data.frame(fread(paste0("GWAS/",fd,"/", sel, "_L.txt.STATUS.glm.logistic"),sep="\t",header=T))

male_df <- data.frame(fread(paste0("GWAS/",fd,"/", sel, "_H.txt.STATUS.glm.logistic"),sep="\t",header=T))

female_df <- female_df[female_df$TEST == "ADD",]
female_df$BETA <- log(female_df$OR)
female_df <- female_df[,c("ID", "A1","BETA","SE")]

male_df <- male_df[male_df$TEST == "ADD",]
male_df$BETA <- log(male_df$OR)
male_df <- male_df[,c("ID", "A1","BETA","SE")]

freq <- freq[freq$ID %in% female_df$ID,]

f1 <- freq[freq$ALT_FREQS < 0.5,]
f1$Minor_allele <- f1$ALT
f1$Major_allele <- f1$REF

f2 <- freq[freq$ALT_FREQS >= 0.5,]
f2$Minor_allele <- f2$REF
f2$Major_allele <- f2$ALT

f <- rbind(f1, f2)

fem<- merge(f, female_df, by = "ID")
fem1 <- fem[fem$Minor_allele == fem$A1,]
fem2 <- fem[fem$Major_allele == fem$A1,]
fem1$BETA <- -fem1$BETA #  polarize by major allele 
femf <- rbind(fem1, fem2)
nrow(femf[femf$BETA > 0,])
nrow(femf[femf$BETA < 0,])

mal<- merge(f, male_df, by = "ID")
mal1 <- mal[mal$Minor_allele == mal$A1,]
mal2 <- mal[mal$Major_allele == mal$A1,]
mal1$BETA <- -mal1$BETA #  polarize by major allele 
malf <- rbind(mal1, mal2)
nrow(malf[malf$BETA > 0,])
nrow(malf[malf$BETA < 0,])

comb <- merge(malf, femf, by = "ID")
comb$BETA.y <- scale(comb$BETA.y, scale =F)
comb$BETA.x <- scale(comb$BETA.x, scale =F)

# create matrix of BETA and SE
# reference covariance matrices: Zhu, et al. Cell Genomics (2023)
conditions <- c("female", "male")
r <- nrow(comb)
BETA <- matrix(c(comb$BETA.y, comb$BETA.x), nrow=r, ncol=2, dimnames=list(c(comb$VAR),conditions))
SE <- matrix(c(comb$SE.y, comb$SE.x), nrow=r, ncol=2, dimnames=list(c(comb$VAR),conditions))

# create mash data object
data = mash_set_data(BETA, SE)
print(head(BETA))
data.temp = mash_set_data(data$Bhat,data$Shat)
data.random = data.temp

# set up canoncial covar matrices and add hypothesis
U.c = cov_canonical(data.random)  
corr = c(1,0.75,0.5,0.25,0,-0.25,-0.5,-0.75,-1)
effect = c(1.5,2,3)
for (c in corr) {
for (e in effect) {
    U.c[[paste('f',c,e,sep="_")]] <- matrix(c(e^2,c*e,c*e,1),2,2)
    U.c[[paste('m',c,e,sep="_")]] <- matrix(c(1,c*e,c*e,e^2),2,2)
}
}
U.c[['equal_-0.25_1']] <- matrix(c(1,-0.25,-0.25,1),2,2)
U.c[['equal_-0.5_1']] <- matrix(c(1,-0.5,-0.5,1),2,2)
U.c[['equal_-0.75_1']] <- matrix(c(1,-0.75,-0.75,1),2,2)
U.c[['equal_-1_1']] <- matrix(c(1,-1,-1,1),2,2)
names(U.c)[1:7] <- c("equal_0_1", "f_0_1", "m_0_1", "equal_1_1", "equal_0.25_1", "equal_0.5_1", "equal_0.75_1")

# fit mash model 
print("start mash")
m = mash(data.random, Ulist= U.c, outputlevel = 1)
print("mashdone")
# mixture model
mixture_prop <- get_estimated_pi(m)
g <- get_fitted_g(m)

mixture <- matrix(names(mixture_prop), ncol=1)
mixture <- cbind(mixture, mixture_prop)
rep <- 1
colnames(mixture) <- cbind(paste0("mix_",0:rep))

df <- mixture
df <- data.frame(df)
df_values <- data.frame(Name = df$mix_0, Mean = df$mix_1, 
                        SE =0)

# split matrice names
df_values <- df_values %>%
separate(Name, c("sex","correlation","magnitude"), sep="[_]", fill="right") %>%
mutate(magnitude = paste0(sex, magnitude))

prepare_df <- function(df) {
df$magnitude <- factor(df$magnitude, levels = c('f1','f3', 'f2', 'f1.5','equal1','m1.5','m2','m3','m1'))
df <- df %>% mutate_at(1, as.numeric) %>%
    arrange(correlation, magnitude)
return(df)
}

# split between null and values
df_ave <- prepare_df(df_values[2:nrow(df_values),c(2,3,4,5)])
df_null <- prepare_df(df_values[1,c(2,3,4,5)])

### SMALL HEATMAP ###
nan_weight <- 1 / (1 - df_null$Mean[1])
df_ave$Mean = df_ave$Mean * nan_weight
df_ave$magnitude <- as.character(df_ave$magnitude) ; df_ave$correlation <- as.numeric(df_ave$correlation)

# group by sex
group_sex <- function(sex){
df_sex <- df_ave %>% filter(substr(magnitude,1,1) == sex) %>%
    group_by(correlation) %>%
    summarise(mean_sum = sum(Mean)) %>%
    as.data.frame()
return(df_sex)
}
for (s in c('f','m','e')) {
assign(s, group_sex(s))
}

# group by correlation
df_small <- data.frame(cbind(e[,1], f[,2], e[,2], m[,2])); colnames(df_small) <- c('corr', 'F', 'E', 'M')
df_small <- data.frame(rbind(
c('perfect', colSums(df_small[df_small$corr == 1, 2:4])),
c('partial,\npositive', colSums(df_small[(df_small$corr > 0) & (df_small$corr < 1), 2:4])),
c('uncorrelated', colSums(df_small[(df_small$corr == 0), 2:4])),
c('negative', colSums(df_small[(df_small$corr < 0), 2:4]))
))
colnames(df_small) <- c('corr', 'L > H', 'L = H', 'L < H')
df_small <- melt(df_small, id.vars=c('corr'))

# table names
df_small$corr <- factor(df_small$corr, levels = c('negative', 'uncorrelated', 'partial,\npositive', 'perfect'))
df_small$variable <- factor(df_small$variable, levels = c('L > H', 'L = H', 'L < H'))
df_small <- df_small %>% mutate_at(3, as.numeric) %>%
arrange(corr, variable)
colnames(df_small) <- c("correlation", "magnitude", "value")

sel1 <- strsplit(sel, "-")[[1]][1]
sel2 <- strsplit(sel, "-")[[1]][2]

des1 <- fielddict$ShortDesc[fielddict$Field == sel1]
des2 <- fielddict$ShortDesc[fielddict$Field == sel2]

df_small$Field <- sel
df_small$F1 <- des1
df_small$F2 <- des2

df_ave$Field <- sel
df_ave$F1 <- des1
df_ave$F2 <- des2


dsmall$correlation <- gsub("[\r\n]", "_", dsmall$correlation)
print("saving files per exp")
dfs <- df_small
dfs$correlation <- gsub("[\r\n]", "_", dfs$correlation)
write.table(dfs, paste0("MASH/perExp/",sel,"_dsmall.txt"), row.names = F, col.names = T, 
        quote= F ,sep="\t")
write.table(df_ave, paste0("MASH/perExp/",sel,"_dbig.txt"), row.names = F, col.names = T, 
        quote= F ,sep="\t")




