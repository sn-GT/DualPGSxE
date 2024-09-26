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

# Prevalence vs PGS curve data

# PGS percentile, prevalence/incidence
PPsave <- PP
table(PPsave$ENVCODE)

# fitted values for low-risk exposure E00
fit1 <- lm(data = PPsave[PPsave$ENVCODE == "00",], Prev ~ poly(PGS, 3, raw = T))
pred1 <- fitted(fit1)
fitted_e1 <- data.frame(x = c(1:100), y = pred1 )

# fitted values for low-risk exposure E11
fit2 <- lm(data = PPsave[PPsave$ENVCODE == "11",], Prev ~ poly(PGS, 3, raw = T))
pred2 <- fitted(fit2)
fitted_e2 <- data.frame(x = c(1:100), y = pred2 )

############################## Area diff ##############################################################################################

# AUC per quartiles 

# Create bins for x values
bins <- cut(fitted_e1$x, breaks = seq(0, 100, by = 10), include.lowest = TRUE) # bins at interval of 10 

selbins <- (list(
  c("[0,10]"),
  c("[0,10]","(10,20]"),
  c("[0,10]","(10,20]", "(20,30]"),
  c("[0,10]","(10,20]", "(20,30]","(30,40]"),
  c("[0,10]","(10,20]", "(20,30]","(30,40]","(40,50]"),
  c("[0,10]","(10,20]", "(20,30]","(30,40]","(40,50]","(50,60]"),
  c("[0,10]","(10,20]", "(20,30]","(30,40]","(40,50]","(50,60]","(60,70]"),
  c("[0,10]","(10,20]", "(20,30]","(30,40]","(40,50]","(50,60]","(60,70]", "(70,80]"),
  c("[0,10]","(10,20]", "(20,30]","(30,40]","(40,50]","(50,60]","(60,70]", "(70,80]" ,"(80,90]"),
  c("[0,10]","(10,20]", "(20,30]","(30,40]","(40,50]","(50,60]","(60,70]", "(70,80]" ,"(80,90]" , "(90,100]")
  ))

# Initialize vectors to store area values for each bin and each group
area_values_e1 <- numeric()
area_values_e2 <- numeric()

# Iterate over bins and calculate area under the curve for each
bin = unlist(bin)

# subset data
subset_data_e1 <- subset(fitted_e1, bins %in% bin)
subset_data_e1 <- na.omit(subset_data_e1)
  
# Check if there are at least two non-NA values for E1
if (sum(!is.na(subset_data_e1$y)) >= 2) {
  # Calculate area under the curve for E1
  area_value_e1 <- integrate(approxfun(subset_data_e1$x, subset_data_e1$y), min(subset_data_e1$x), max(subset_data_e1$x))$value
  area_values_e1 <- c(area_values_e1, area_value_e1)
} else {
  area_values_e1 <- c(area_values_e1, NA)
}
  
# Extract data for the specific bin for E2
#subset_data_e2 <- subset(fitted_e2, bins == bin)
subset_data_e2 <- subset(fitted_e2, bins %in% bin)
subset_data_e2 <- na.omit(subset_data_e2)
  
# Check if there are at least two non-NA values for E2
if (sum(!is.na(subset_data_e2$y)) >= 2) {
  # Calculate area under the curve for E2
  area_value_e2 <- integrate(approxfun(subset_data_e2$x, subset_data_e2$y), min(subset_data_e2$x), max(subset_data_e2$x))$value
  area_values_e2 <- c(area_values_e2, area_value_e2)
} else {
  area_values_e2 <- c(area_values_e2, NA)
}
}

# Create data frames with bin labels and corresponding area values for each group
area_data_e1 <- data.frame(bin = levels(bins), area_e1 = area_values_e1)
area_data_e2 <- data.frame(bin = levels(bins), area_e2 = area_values_e2)

# Display the area values for each bin and each group
print("Area under curve for E1:")
print(area_data_e1)

print("Area under curve for E2:")
print(area_data_e2)

area_data_e1$area_e1 <- area_data_e1$area_e1/100
area_data_e2$area_e2 <- area_data_e2$area_e2/100

areadf <- merge(area_data_e1, area_data_e2, by = "bin")

# Proportion needed to benefit (100/AUC difference between low and high-risk curves)
areadf$areadiff <- abs(areadf$area_e1 - areadf$area_e2)
areadf$Prop <- 100/areadf$areadiff

###################################################################################################################################################################

# Plotting area diff and proportion benefit  (cumulative )

ord <- c("[0,10]" ,  "(10,20]" , "(20,30]" , "(30,40]",  "(40,50]" , "(50,60]" , "(60,70]" , "(70,80]",  "(80,90]",  "(90,100]")
areadf$bin <- factor(areadf$bin, levels = ord)

areadf <- areadf[order(factor(areadf$bin, levels = ord)), ]
areadf$Scale_num <- c(1:10)
areadf$Scale_den <- c(10)
areadf$scaled <- (areadf$Prop *areadf$Scale_num)/areadf$Scale_den

areadf$disc <- round(areadf$scaled, 0)
uq <- areadf[!duplicated(areadf$disc),]
areadf$disc2 <- areadf$disc
areadf$disc2[!areadf$bin %in%  uq$bin] <- ""

areadf$field1 <- env_selected[1]
areadf$field2 <- env_selected[2]
areadf$env1 <- field_dict$Description[field_dict$Field == env_selected[1]]
areadf$env2 <- field_dict$Description[field_dict$Field == env_selected[2]]
areadf$des1 <- field_dict$ShortDesc[field_dict$Field == env_selected[1]]
areadf$des2 <- field_dict$ShortDesc[field_dict$Field == env_selected[2]]

nnt <- rbind(nnt, areadf)
###########################################################################################################################################################

# Add PNB bar below prev-PGS curves 

gg_diff_legend <- ggplot(areadf, aes(x = bin, y = 0.5, fill = scaled)) +
  geom_tile() +
  scale_fill_gradient(low = "lightblue", high = "maroon", name = "Absolute Difference") +
  theme_minimal() +
  theme(
    legend.position = c(0.5, 0.85),  # Adjust the position of the legend
    legend.direction = "horizontal",  # Display legend items horizontally
    legend.key.size = unit(1, "cm"),  # Adjust the size of legend items
    axis.title = element_blank(),
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(size = 11),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  scale_x_discrete(breaks =levels(areadf$bin), labels = areadf$disc2) +  # Use scale_x_discrete for discrete x values
  coord_cartesian(ylim = c(0, 1))  # Set y-axis limits to avoid cutting off the plot

gg_diff_legend

# save figure
png(paste0("//Plots_PNB/", env_selected[1],"-",env_selected[2],".png"), units="in", width=5, height=5.5, res=200)
p <- plot_grid(p1, gg_diff_legend + theme(legend.position = "None"), ncol = 1, align = 'v', rel_heights = c(2, 0.2))
print({p})
dev.off()
