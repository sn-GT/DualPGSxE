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

# pairwise combinations of exposures
checking 2
# select pairwise combinations