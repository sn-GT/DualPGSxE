library(data.table)
library(ggplot2)
library(cowplot)
library(scales)
library(plyr)
require(gridExtra)
require(grid)
library(ggpmisc)
library(gghighlight)
library(Rmisc)

# Load functions 
source(paste0("/functions/prevperc_func.R"))
source(paste0( "/functions/additive_CompEnv_func_PrevLH.R"))
source(paste0("/functions/interaction_CompEnv_func.R"))
source(paste0("/functions/delta_CompEnv_rldev.R"))


trait <- "CAD"
xl <- "CAD"

# Polygenic Score (PGS)
pgs <- read.table("PGS_CAD.sscore", header = T)
pgs <- pgs[,c("IID", "SCORE1_SUM")]
names(pgs)[2] <- "PGS"

# Case control status for disease
# Subset for white british
cc <- read.table("CAD_CaseCont_WB.txt", header = T)
names(cc) <- c("IID", "CC")

# merge pgs and case-control
combined_df <- merge(pgs, cc, by = "IID")

# selected field for exposures in the UKB 
# field id, exposure name, category, sub-category
field_dict <- read.table("data/field_dict.txt", header = T)
fields <- field_dict$FieldID

# pairwise combinations of exposures
exposure_pairs <- data.frame(t(combn(cont,2)))
# columns: FieldID1, FieldID2
names(exposure_pairs) <- c("Exposure1_fieldID", "Exposure2_fieldID")

# for example selecting first exposure pair
select_exposure <- exposure_pairs[1,]

# Get phenotype data for exposure pairs
# columns IID, phenotype, GROUP
phe1 <- read.table(paste0("/UKB_categories/", select_exposure$Exposure1_fieldID, ".txt"), header = T, sep = "\t")
phe2 <- read.table(paste0("/UKB_categories/", select_exposure$Exposure2_fieldID, ".txt"), header = T, sep = "\t")
pheno_sub <- merge(phe1, phe2, by = "IID") # check nrow(pheno_sub) != 0

pheno_sub$CompGROUP <- paste0(pheno_sub$GROUP.x, "|", pheno_sub$GROUP.y)
table(pheno_sub$CompGROUP) # check length(table(pheno_sub$CompGROUP)) == 4

combined_sub <- combined_df[,c("IID", "PGS", "CC")]

# merge 
combined_sub <- merge(combined_sub, pheno_sub, by = "IID")

# check prevalence in each group
cp <- summarySE(combined_sub, measurevar = "CC", groupvars = c("CompGROUP"))
cp <- cp[order(cp$CC, decreasing =T),]

pheno_sub$ENVCODE[pheno_sub$CompGROUP == cp$CompGROUP[1]] <- "11"
pheno_sub$ENVCODE[pheno_sub$CompGROUP == cp$CompGROUP[2]] <- "10"
pheno_sub$ENVCODE[pheno_sub$CompGROUP == cp$CompGROUP[3]] <- "01"
pheno_sub$ENVCODE[pheno_sub$CompGROUP == cp$CompGROUP[4]] <- "00"

pheno_sub$CompGROUP <- paste0(pheno_sub$CompGROUP, "\n", pheno_sub$ENVCODE)

pheno_majorgroup <- pheno_sub[,c("IID", "CompGROUP", "ENVCODE")]

m <- merge(pheno_majorgroup, combined_sub, by = "IID")
names(m) <- c("FID", "GROUP","ENVCODE", "PGS", "CC")

# Get group names of combined DF
l <- names(table(m$GROUP))

# get overall prevalence of combined DF
overallprev <- table(m$CC)[2]/nrow(m)

# compute prevalence vs percentile PGS
PP <- prevperc(m)

# order by prevalence and assign color
prev_groups <- summarySE(PP, measurevar = "Prev", groupvars = c("GROUP"))
prev_groups$ENVCODE <- gsub(".*\n", "", prev_groups$GROUP) 
prev_groups <- prev_groups[order(prev_groups$Prev, decreasing = T),]
prev_groups$COL <- c("#7B6079",  "#DE8971","#A7D0CD", "#e6c78e")
col <- as.character(prev_groups$COL)
names(col) <- as.character(prev_groups$GROUP)

PPsave <- PP
          
pattern <- c("11", "10", "01", "00")  # Desired patterns in the desired order

# Function to extract the pattern
extract_pattern <- function(x) {
for (p in pattern) {
  if (grepl(p, x))
    return(p)
}
return("")  # Handle cases where no match is found (though it shouldn't happen in this case)
}

# Apply the function to extract the pattern
PPsave$pattern <- sapply(PPsave$GROUP, extract_pattern)
PPsave$GROUP <- factor(PPsave$GROUP, levels = unique(PPsave$GROUP[order(match(PPsave$pattern, pattern))]))

# PLot prevalence vs percentile curve
p1 <- ggplot(PPsave, aes(x=PGS, y=Prev, col = GROUP, group = GROUP)) + 
geom_point( size = 1.2) + theme_bw() + 
xlab(paste0("Percentile of PGS-", xl)) + 
ylab(paste0("Incidence of ",trait, "(%)")) + 
ggtitle(paste0(field_dict$Description[field_dict$Field == env_selected[1]], " -\n", field_dict$Description[field_dict$Field == env_selected[2]])) + 
geom_smooth(method = "lm", formula = y ~ poly(x, 3, raw = T), size = 1.2, se = F) + theme(legend.position = "top") +  
theme_classic() + 
theme(legend.position = "top", legend.box = "vertical", legend.margin=margin(),
      legend.text = element_text(size=12, color="black"),
      title = element_text(size = 13), 
      axis.text.x = element_text( size = 12, color = "black"), 
      axis.title.x  = element_text(size = 12, color = "black"),
      axis.title.y  = element_text(size = 12, color = "black"),
      axis.text.y  = element_text(size = 12, color = "black")) +  
scale_color_manual(name = "", values = col) + 
guides(color=guide_legend(nrow=2,byrow=TRUE)) + 
theme(
  panel.background = element_rect(fill='transparent'),
  plot.background = element_rect(fill='transparent', color=NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.background = element_rect(fill='transparent'),
  legend.box.background = element_rect(fill='transparent', color=NA)
) + ylim(-0.5,38)
p1

###############################################################################################

# Additive model - expected curves 

PP_full <- add_comp(PP, overallprev)

# Compute R2 of the model comparing expected to observed curves 
r2 <- 1 - (var(PP_full$Prev - PP_full$ExpectedPrev_nostoch)/ var(PP_full$Prev))
n=400
k=3
adjr2 <- 1 - ((1 -r2)*(n-1)/(n-k-1))
PP_full$AdditiveR2_prev <- adjr2
PP_full$field1 <- paste0(env_selected[1])
PP_full$field2 <- paste0(env_selected[2])
PP_full$Envdesc <- paste0(field_dict$Description[field_dict$Field == env_selected[1]], " - ", field_dict$Description[field_dict$Field == env_selected[2]])

# Get average expected prevalence + std error from 50 iterations
PP_summ1 <- summarySE(PP_full, measurevar=c("ExpectedPrev"), groupvars=c("PGS", "GROUP", "ENVCODE"))
PP_summ1$GROUP[PP_summ1$GROUP == l[1]] <- paste0(l[1], "_expected")
PP_summ1$GROUP[PP_summ1$GROUP == l[2]] <- paste0(l[2], "_expected")
PP_summ1$GROUP[PP_summ1$GROUP == l[3]] <- paste0(l[3], "_expected")
PP_summ1$GROUP[PP_summ1$GROUP == l[4]] <- paste0(l[4], "_expected")
names(PP_summ1)[5] <- "Prev"
PP_summ1$ENVCODE <- paste0(PP_summ1$ENVCODE, "_expected")

PP_summ2 <- summarySE(PP_full, measurevar="Prev", groupvars=c("PGS", "GROUP", "ENVCODE"))
PP_summ <- rbind(PP_summ1, PP_summ2)

col2 <- c(col, rep("grey60", 4))
names(col2) <- c(names(col), paste0(names(col), "_expected"))

lt <- c(rep("solid",4), rep("longdash", 4))
names(lt) <- c(names(col), paste0(names(col), "_expected"))

# Prevalence vs percentile curve with overlaid grey curves for expected risk under null model of additive expectation
p2 <- ggplot(PP_summ, aes(x=PGS, y=Prev, group = GROUP, color = GROUP, linetype=GROUP)) + 
xlab(paste0("Percentile of PGS-", xl)) + 
ylab(paste0("Incidence of ",trait)) + 
theme_classic() + 
ggtitle(paste0("Additive R2 = ",round(PP_full$AdditiveR2[1], 2)),
        subtitle = bquote("ui ~ " * bar(PGS) * " + E1 + E2") ) +
geom_smooth(method = "lm", formula = y ~ poly(x, 3, raw = T), size = 1.2, se = F)+
theme(legend.position = "None", legend.box = "vertical", legend.margin=margin(),
      title = element_text(size = 11), 
      axis.text.x = element_text( size = 12, color = "black"), 
      axis.title.x  = element_text(size = 12, color = "black"),
      axis.title.y  = element_text(size = 12, color = "black"),
      axis.text.y  = element_text(size = 12, color = "black")) +  
scale_color_manual(name = "", values = col2) + 
scale_linetype_manual(name="",values=lt)  + 
theme(
  panel.background = element_rect(fill='transparent'),
  plot.background = element_rect(fill='transparent', color=NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.background = element_rect(fill='transparent'),
  legend.box.background = element_rect(fill='transparent', color=NA)
)
p2
###############################################################################################

# Delta : deviations in disease risk at the extremes of PGS for observed vs expected curves under null model

codes <- prev_groups$ENVCODE
envcomb <- list(c(codes[4], codes[3]), c(codes[4], codes[2]),  c(codes[4], codes[1]), 
              c(codes[3], codes[2]), c(codes[3], codes[1]), 
              c(codes[2], codes[1]))
          
# function to compute delta           
delfull <- delta_comp(envcomb, PP_full)

adddelfull <- delfull
deladd <- rbind(deladd, adddelfull)

# observed: delta at top 2 percentile, PGS > or < 2sd
delobs2 <- summarySE(delfull, measurevar = "deltaObs2", groupvars = c("Delta", "ENV"))
delobs2$Model <- "Observed_top2"
names(delobs2)[4] <- "DeltaVal"

# observed: delta at top 5 percentile, PGS > or < 5sd
delobs5 <- summarySE(delfull, measurevar = "deltaObs5", groupvars = c("Delta", "ENV"))
delobs5$Model <- "Observed_top5"
names(delobs5)[4] <- "DeltaVal"

# expected: delta at top 2 percentile, PGS > or < 2sd
delexp2 <- summarySE(delfull, measurevar = "deltaExp2", groupvars = c("Delta", "ENV")) 
delexp2$Model <- "Exp_Additive_top2"
names(delexp2)[4] <- "DeltaVal"

# expected: delta at top 5 percentile, PGS > or < 5sd
delexp5 <- summarySE(delfull, measurevar = "deltaExp5", groupvars = c("Delta", "ENV")) 
delexp5$Model <- "Exp_Additive_top5"
names(delexp5)[4] <- "DeltaVal"

delexp2_nostoc <- summarySE(delfull, measurevar = "deltaExp2_nostoc", groupvars = c("Delta", "ENV")) 
delexp2_nostoc$Model <- "Exp_Additive_top2_nostoc"
names(delexp2_nostoc)[4] <- "DeltaVal"

###############################################################################################

# Interaction model - expected curves 

  PP <- PPsave
  l <- names(table(PP$GROUP))

  # Expected risk under the interaction model
  PP_full <- int_comp(PP, overallprev)

  # Compute R2 of observed vs expected risk under interaction model 
  r2 <- 1 - (var(PP_full$Prev - PP_full$ExpectedPrev_nostoch)/ var(PP_full$Prev))
  n=400
  k=3
  adjr2 <- 1 - ((1 -r2)*(n-1)/(n-k-1))
  PP_full$InteractionR2_prev <- adjr2
  PP_full$field1 <- paste0(env_selected[1])
  PP_full$field2 <- paste0(env_selected[2])
  PP_full$Envdesc <- paste0(field_dict$Description[field_dict$Field == env_selected[1]], " - ", field_dict$Description[field_dict$Field == env_selected[2]])
  
  interactionPP <- rbind(interactionPP, PP_full[1,])
  
  # Get average expected prevalence + std error from 50 iterations
  PP_summ1 <- summarySE(PP_full, measurevar="ExpectedPrev", groupvars=c("PGS", "GROUP", "ENVCODE"))
  PP_summ1$GROUP <- paste0(PP_summ1$GROUP, "_expected")
  PP_summ1$ENVCODE <- paste0(PP_summ1$ENVCODE, "_expected")
  names(PP_summ1)[5] <- "Prev"
  
  PP_summ2 <- summarySE(PP_full, measurevar="Prev", groupvars=c("PGS", "GROUP", "ENVCODE"))
  
  PP_summ <- rbind(PP_summ1, PP_summ2)
  interactionPPsave <- PP_summ
  
  # Prevalence vs percentile curve with overlaid grey curves for expected risk under null model of additive expectation
  p3  <- ggplot(PP_summ, aes(x=PGS, y=Prev, group = GROUP, color = GROUP, linetype=GROUP)) + 
    xlab(paste0("Percentile of PGS-", xl)) + 
    ylab(paste0("Incidence of ",trait)) + 
    theme_classic() + 
    ggtitle(paste0("Interaction R2 = ",round(PP_full$InteractionR2[1], 2)),
            subtitle = bquote("ui ~ " * bar(PRS) * " + E1 + E2 +" * bar(PRS) * "xE1 + " * bar(PRS) * "xE2 + " * bar(PRS) * "xE1xE2") ) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3, raw = T), size = 1.2, se = F)+
    theme(legend.position = "None", legend.box = "vertical", legend.margin=margin(),
          title = element_text(size = 11), 
          axis.text.x = element_text( size = 12, color = "black"), 
          axis.title.x  = element_text(size = 12, color = "black"),
          axis.title.y  = element_text(size = 12, color = "black"),
          axis.text.y  = element_text(size = 12, color = "black")) +  
    scale_color_manual(name = "", values = col2) + 
    scale_linetype_manual(name="",values=lt)  + #ylim(0, 0.7) +  
    theme(
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill='transparent'),
      legend.box.background = element_rect(fill='transparent', color=NA)
    )
  p3


###############################################################################################

# compute delta of observed vs expected risk under the interaction model

delfull <- delta_comp(envcomb, PP_full)
intdelfull <- delfull
delint <- rbind(delint, intdelfull)

delint2 <- summarySE(delfull, measurevar = "deltaExp2", groupvars = c("Delta", "ENV")) 
delint2$Model <- "Exp_Interaction_top2"
names(delint2)[4] <- "DeltaVal"

delint5 <- summarySE(delfull, measurevar = "deltaExp5", groupvars = c("Delta", "ENV")) 
delint5$Model <- "Exp_Interaction_top5"
names(delint5)[4] <- "DeltaVal"

delint2_nostoc <- summarySE(delfull, measurevar = "deltaExp2_nostoc", groupvars = c("Delta", "ENV")) 
delint2_nostoc$Model <- "Exp_Interaction_top2_nostoc"
names(delint2_nostoc)[4] <- "DeltaVal"

###############################################################################################

