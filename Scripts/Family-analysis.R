
# Part 4 (family analysis)

#======================================================================================================

#Load required packages
library(Hmisc)
library(tidyverse)
library(reshape2)
library(here)
library(xlsx)
library(dplyr)
library(ggformula)
library(nlme)

#======================================================================================================

# Get number of species in each family across the dataset
# It needs the original community data with two first extra columns
# First extra column: family; second extra column: genus (in case of analysing genera data)

family <- openxlsx::read.xlsx(here("Data", "data.family.xlsx"))

# Column family as factor
family$family <- factor(family$family)

# Function to get a family X site matrix
# Kindly provided by dr. Aurèle Toussaint
family_by_site <- matrix(0, nc = ncol(family) - 3, nr = length(levels(family[, 1])),
                         dimnames = list(levels(family[, 1]), colnames(family)[ -c(1 : 3)]))

for (site in 4: ncol(family)){
  fam_site <- family[, c(1, site)]
  family_by_site[names(table(fam_site[which(fam_site[, 2] == 1), 1])), site - 3] <- table(fam_site[which(fam_site[, 2] == 1), 1])
}


#Get a transposed matrix
site.by.family <- as.data.frame(t(family_by_site))

#======================================================================================================

# Get number of species per family
number.of.species.per.family <- family %>% 
  group_by(family) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n))

#======================================================================================================
# Load dataset
data <- openxlsx::read.xlsx(here("Data", "data.full.xlsx"))

# Rename eco.region column from the dataset and transform it into factor
data$eco.region <- as.factor(data$eco.region)

levels(data$eco.region)[match("Northeastern Atlantic Shorelands",levels(data$eco.region))] <- "Northeastern restingas"
levels(data$eco.region)[match("Central Atlantic Shorelands",levels(data$eco.region))] <- "Central restingas"
levels(data$eco.region)[match("Southern Atlantic Shorelands",levels(data$eco.region))] <- "Southern restingas"
data$eco.region <- factor(data$eco.region, levels = c("Northeastern restingas", "Central restingas", "Southern restingas"))


#=====================================================================================================
# Get 10 most abundant families in the site.by.family dataset and plot them

myrtaceae<- site.by.family$myrtaceae
fabaceae <- site.by.family$fabaceae
rubiaceae <- site.by.family$rubiaceae
lauraceae <- site.by.family$lauraceae
melastomataceae <- site.by.family$melastomataceae
sapotaceae <- site.by.family$sapotaceae
annonaceae <- site.by.family$annonaceae
euphorbiaceae <- site.by.family$euphorbiaceae
moraceae <- site.by.family$moraceae
solanaceae <- site.by.family$solanaceae

# Put it into a data frame
most.abund.families <- data.frame(myrtaceae, fabaceae, rubiaceae, lauraceae,
                                  melastomataceae, sapotaceae, annonaceae, euphorbiaceae, moraceae, solanaceae)

# Rownames - including sites and latitude
rownames(most.abund.families) <- rownames(site.by.family)
most.abund.families$lat <- data$lat
most.abund.families$lon <- data$lon
most.abund.families$richness <- data$species
most.abund.families$site <- rownames(site.by.family)
most.abund.families$eco.region <- data$eco.region
colnames(most.abund.families) <-   capitalize(colnames(most.abund.families))

# Transform the dataset into a format to plot everything together
All.families <- melt(most.abund.families, variables.name = c("Myrtaceae", "Fabaceae", "Rubiaceae",   "Lauraceae", "Melastomataceae", 
                                                             "Sapotaceae", "Annonaceae", "Euphorbiaceae", "Moraceae", "Solanaceae"),
                     id = c("Site", "Lat", "Lon", "Eco.region", "Richness"))

#=====================================================================================================

# Plot the result
tiff("most.abundant.families.tiff", units="in", width=14, height=7, res=600)

# Species per family along the latitudinal gradient
family.plot <- ggplot(All.families, aes(abs(Lat), value))+
  geom_point(aes(color = Eco.region))+
  geom_lm(aes(abs(Lat), value), formula = y~poly(x, 2), color = "dark blue", lwd = 1.2)+
  facet_wrap(~variable, nrow = 2, ncol = 5)+
  theme_bw()+
  scale_y_continuous(name = "Species richness / family")+
  scale_x_continuous(name = "Latitude")+
  theme(axis.title.y = element_text(size = 17))+
  theme(axis.title.x = element_text(size = 17))+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.text.x = element_text(size=12))+
  theme(strip.text = element_text(size = 13))+
  scale_color_viridis_d()+
  theme(panel.border = element_rect(fill = NA, colour = "black"))

# Close the plot window
dev.off()

# NOTE: As we wrote in part 3, we plotted the figure 3 with 95% individually
# for each of the ten families. Then we grouped them into one big panel.
# We will provide the code for one family, and you can use it the other nine.
# If you know a better (more elegant/quick) way to do it, let me know. It would
# be nice to see it.

# GLS model and figure of the
# latitudinal gradient in species richness of the richnest
# angiosperm family in restingas - Myrtaceae


#=====================================================================================================

# Fit the model
mod.myrtaceae <- gls(Myrtaceae ~ poly(abs(Lat), 2),
                     correlation = corExp(form = ~Lon + Lat), data = most.abund.families)

# Get predicted values
fit <- predict(mod.myrtaceae)

# Variace-covariance matrix
V <- vcov(mod.myrtaceae)

# Model matrix corresponding to the new data
X <- model.matrix( ~poly(abs(lat), 2),data = data)

# Standard errors
se.fit <- sqrt(diag(X %*% V %*% t(X)))

# Data frame for the confidence interval
pred <- with(most.abund.families, data.frame(Lat,
                                                  Myrtaceae = fit, lwr = fit-1.96*se.fit, upr = fit+1.96*se.fit))


# Plot it
myrtaceae.plot <- ggplot(most.abund.families, aes(abs(Lat), Myrtaceae))+
  geom_point(aes(color = data$eco.region))+
  geom_line(data=pred, color = "dark blue", lwd = 1.2, alpha = 0.5)+
  geom_ribbon(data=pred,aes(ymin=lwr,ymax=upr), alpha=0.2)+
  scale_x_continuous(name = NULL)+
  scale_y_continuous(name = NULL)+
  xlab(NULL)+
  ylab(NULL)+
  theme(axis.title.y = element_text(size = 17))+
  theme(axis.title.x = element_text(size = 17))+
  scale_color_viridis_d()+
  theme_bw()+
  labs(color = "")+
  guides(colour = guide_legend(override.aes = list(size=4)))+
  ggtitle(label = "Myrtaceae")

# END
