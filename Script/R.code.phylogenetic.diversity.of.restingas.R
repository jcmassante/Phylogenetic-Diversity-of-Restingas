# =======================================================================================================
# R code used in the following publication:
# =======================================================================================================
# Title:   ENVIRONMENT AND EVOLUTIONARY HISTORY DEPICT PHYLOGENETIC ALPHA AND BETA DIVERSITY
# IN THE ATLANTIC COASTAL WHITE-SAND WOODLANDS

# Authors: Jhonny C Massante & Pille Gerhold
#======================================================================================

# Standardisation of species names
# Phylogenetic tree
# Indexes of phylogenetic alpha and beta diversity
# Statistical analysis of main and supplementary results
# Figures of main and supplementary texts
#======================================================================================

# All analyses were carried out using R Project. Users will need to 
# create their own R Project or working directory.

#======================================================================================

# Load required packages
library(Taxonstand)
library(V.PhyloMaker)
library(ape)
library(picante)
library(phytools)
library(betapart)
library(openxlsx)
library(mapdata)
library(spdep)
library(geosphere)
library(tidyverse)
library(here)
library(corrplot)
library(nlme)
library(MuMIn)
library(MASS)
library(ecodist)
library(car)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggpubr)
library(raster)
library(gridExtra)
library(grid)
library(ggsn)


#======================================================================================
# Standardisation of species names
# Input: data frame with original (from community data) species names according to the Taxonstand package requirements
# Output: data frame with standardised species names
# NOTE: it may be necessary to do manual standardisation of some specific taxa

to.standardise <- read.xlsx(here("Data", "to.standardise.xlsx"))# Input
names <- TPL(to.standardise$Full.name)# Output

# Save the output
write.xlsx(names, "splist.xlsx")


#======================================================================================

# Phylogenetic tree
# Input: data drame with species list (standardised names) according to V.Phylomaker package requirements
# Output: a time calibrated phylogenetic tree

# Read required species list
splist <- read.xlsx(here("Data", "splist.xlsx"))# Input

# Building tree using phylo.maker function from V.Phylomaker package
Phylogenetic.tree <- phylo.maker(sp.list = splist, tree = GBOTB.extended, nodes = nodes.info.1, scenarios = "S3")# Output

# Exclude tree fern species (there was a gymnosperm species ("podocarpus sellowii") which has been excluded manually)
tree.fern <- c("Cyathea_microdonta" , "Cyathea_leucofolis", "Cyathea_delgadii", "Cyathea_corcovadensis", "Cyathea_atrovirens")
tree <- drop.tip(tree, tree.fern)


#======================================================================================

# Calculation of response variables (phylogenetic alpha and beta diversity)


# Calculate phylogenetic alpha diversity as the standardised effect size of the
# mean pairwise distance (MPD) and mean nearest taxon distance (MNTD) using the Picante package
# Input: matrix of community data (species x sites) and a phylogenetic tree


# Phylogenetic tree
tree <- read.tree(here("Data", "Phylogenetic.tree"))# Input

# Community data
comm <- read.xlsx(here("Data", "comm.xlsx"), rowNames = T)# Input
colnames(comm) <- gsub("\\.", "-", colnames(comm))# Replace symbols to prevent error

# Prune the phylogenetic tree
tree.pruned<- prune.sample(comm, tree)

# Sort the columns of community data to be in the same order as the tip labels of the phylogeny
comm<- comm[, tree.pruned$tip.label]

# Create species distance matrix from the pruned tree
dist<-cophenetic(tree)


### Phylogenetic alpha diversity ###

# Calculate standardised effect size of mean pairwise distance (SESmpd)
sesmpd<- ses.mpd(as.matrix(comm), (dist), null.model = "taxa.labels", runs = 999)# Output
#write.table(sesmpd, "sesmpd.csv", sep=",")# To include in the main dataset

# Calculate standardised effect size of mean nearest taxon distance (SESmntd)
sesmntd <- ses.mntd(as.matrix(comm), (dist), null.model = "taxa.labels", runs = 999)# Output
#write.table(sesmntd, "sesmntd.csv", sep=",")# To include in the main dataset

#======================================================================================

### Phylogenetic beta diversity ###

# Phylogenetic beta diversity metrics from Picante package.

# Comdist calculates the expected phylogenetic distance
# separating two individuals or taxa drawn randomly from different communities.
# Comdistnt calculates the average phylogenetic distance to the most similar 
# taxon or individual in the other community for taxa or individuals in two communities.

# Input: distance matrix calculated from the phylogenetic tree (line 86)

# Deep phylogenetic levels
phylobeta.dpw <- comdist(comm, dist)
save(phylobeta.dpw, file = "phylobeta.dpw.RData")# Save for further analysis

# Shallow phylogenetic levels
phylobeta.dnn <- comdistnt(comm, dist)
save(phylobeta.dnn, file = "phylobeta.dnn.RData")# Save for further analysis


# Phylogenetic beta diversity from Betapart package
# phylo.beta.pair calculates dissimilarity matrix accounting for phylogenetic beta diversity,
# measured as Sorensen derived pairwise phylogenetic dissimilarity.

phylo.restinga <- phylo.beta.pair(comm, tree)# Total phylobeta (turnover + nestedness)
simpson.true.phylo.beta <- phylo.restinga$phylo.beta.sim# Only turnover (Simpson which is unaffected by species richness)
save(simpson.true.phylo.beta, file = "simpson.true.phylo.beta.RData")# Save for further analysis


#======================================================================================

# Brazilian map with restinga assemblages (main figure 1)

# The map included in this study was created following a tutorial available at
# https://rstudio-pubs-static.s3.amazonaws.com/176768_ec7fb4801e3a4772886d61e65885fbdd.html

# Read shapefile
br <- rgdal::readOGR(dsn = "", layer = "Brazil_biomes")

# Check projection
proj4string(br)

# Check names
names(br)
head(br$name)

# Select atlantic biome                     
atlantic <- br[br@data$name=="Mata AtlÃ¢ntica",]
MA <- atlantic

# Plot atlantic biome in ggplot2
brazil.map <- borders("worldHires", regions = "Brazil", fill = "white", colour = "black")
bra.map <- map_data("worldHires","Brazil")
map.1 <- ggplot() + brazil.map + coord_equal()
MAf <- fortify(MA)#Fortify to bw used in ggplot2

# Rename eco.region column from the dataset and transform it into factor
data$eco.region <- as.factor(data$eco.region)

levels(data$eco.region)[match("Northeastern Atlantic Shorelands",levels(data$eco.region))] <- "Northeastern restingas"
levels(data$eco.region)[match("Central Atlantic Shorelands",levels(data$eco.region))] <- "Central restingas"
levels(data$eco.region)[match("Southern Atlantic Shorelands",levels(data$eco.region))] <- "Southern restingas"
data$eco.region <- factor(data$eco.region, levels = c("Northeastern restingas", "Central restingas", "Southern restingas"))

# Fortify data.frame with geographic coordinates of sites (taken from the dataset)
coord <- data[, 8:9]
coord <- fortify(coord)

# Re-order the eco.region factor to be from Northeastern to Southern and not alphabetically
coord$eco.region <- factor(data$eco.region, levels = c("Northeastern restingas", 
                                                       "Central restingas",
                                                       "Southern restingas"))
# Plot brazil map with atlantic biome
map.2 <- map.1 +
  geom_polygon(data = MAf, aes(x = long, y = lat, group = group), fill = "grey", alpha = 0.5)+
  theme_bw()
#======================================================================================

# Get South American map

# Global map
world <- map_data("world")

# A vector with all countries of South America
countrylist <- read.table(here("Data", "countrylist.txt"))

# Get coordinates of countries
lat <- map_data("world",region=countrylist$V1)

# Remove some islands for visual effects
idx <- which(lat$long<=-75 & lat$lat<=-20)

# Plot South American map
g1 <- ggplot() + geom_polygon(data = lat[-idx,], aes(x=long, y = lat, group = group),fill='white', colour="black") +
  coord_equal()+
  theme_bw()+
  geom_polygon(data = MAf, aes(x = long, y = lat, group = group), fill = "dark grey", alpha = 0.8)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y =element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())

# Plot Brazilian map with the Atlantic forest domain and assemblages used in this study

# Plot the map
g2 <- map.2+
  geom_point(data = coord, aes(lon, lat, color = eco.region), size = 1.1, alpha = 0.8)+
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  labs(color = "")+
  scale_colour_viridis_d()+
  ggsn::scalebar(bra.map, dist = 500, dist_unit = "km", 
                 transform = T, st.size=2.5, height=0.02, model = 'WGS84')+
  north(bra.map)+
  theme(axis.text.x = element_text(size = 17))+
  theme(axis.text.y = element_text(size = 17))+
  theme(legend.text = element_text(size = 17))

# Plot both maps together (South American map as an inset map) and export them as a tiff file
tiff(file="mapa.tiff", units = "in", w= 12, h=9, res=600)
grid.newpage()
v1<-viewport(width = 1, height = 1.2, x = 0.5, y = 0.5)# Plot area for the main map
print(g2,vp=v1) 

v2<-viewport(width = 0.28, height = 0.3, x = 0.15, y = 0.32)# Plot area for the inset map
print(g1,vp=v2)

# Close the plot window
dev.off()



#======================================================================================
#Preparation of predictors for analysis of phylogenetic alpha diversity

#Load dataset
data<-openxlsx::read.xlsx(here("Data", "data.full.xlsx"))

# Rename soil.fertility and transform it into factor
levels(data$soil.fertility)[match("Dystrophic (10-25% TBS)",levels(data$soil.fertility))] <- "Dystrophic"
levels(data$soil.fertility)[match("Hypo-Dystrophyc (0-10% TBS)",levels(data$soil.fertility))] <- "Hypo"
levels(data$soil.fertility)[match("Mesotrophic (25-50% TBS)",levels(data$soil.fertility))] <- "Mesotrophic"
data$soil.fertility <- factor(data$soil.fertility, levels = c("Hypo", "Dystrophic", "Mesotrophic"))

# Rename soil.salinity and transform it into factor
levels(data$soil.salinity)[match("Low salinity (ECe ca. 0.3-0.4 ds/m)",levels(data$soil.salinity))] <- "Low.salinity"
levels(data$soil.salinity)[match("Slightly saline (ECe ca. 0.1-0.2 ds/m)",levels(data$soil.salinity))] <- "Slightly.salinity"
data$soil.salinity <- factor(data$soil.salinity, levels = c("Low.salinity", "Slightly.salinity"), ordered = F)


# Variables to check the correlation
vars.mpd <- c("lat", "annual.mean.temp", "temp.seasonality", "annual.precip", "precip.seasonality", "max.temp.warmest.month",
              "precip.driest.month", "clim.stab.temp", "clim.stab.precip")

# Perform pairwise correlation
data.mpd <- data[, vars.mpd]
corr.mpd <- cor(data.mpd)

# Prepare to save the figure (Appendix S1, Fig.S1.1)
tiff("correlation.tiff", units="in", width=10, height=10, res=600)

# Plot the figure
corrplot(corr.mpd, type = "upper", diag = F, method = "number",
         number.cex = 1.5, number.font = 2, tl.cex = 1.5, tl.col = "black")

# Close the plot window 
dev.off()

#Excluded variables strongly correlated with latitude and limiting factors(r > 0.60) . It included
#annual mean temperature, temp.seasonality, precip seasonality, annual precipitation. 
vars.mpd <- c("lat", "max.temp.warmest.month", 
              "precip.driest.month", "clim.stab.temp", "clim.stab.precip", "sesmpd.qian")


# Scale the variables, except latitude
data.mpd <- data[, vars.mpd]
data.mpd <- scale(data.mpd[, -1])

# Put latitude back and include categorical variables (soil fertility and soil salinity)
data.mpd <- as.data.frame(data.mpd)
data.mpd$lat <- abs(data$lat)
data.mpd$lon <- data$lon
data.mpd$soil.fertility <- as.factor(data$soil.fertility)
data.mpd$soil.salinity <- data$soil.salinity


#======================================================================================
# Repeat the procedure for SESmntd

vars.mntd <- c("lat", "max.temp.warmest.month", 
               "precip.driest.month", "clim.stab.temp", "clim.stab.precip", "sesmntd.qian")

# Get the variables from the data frame and scale them
data.mntd <- data[, vars.mntd]
data.mntd <- scale(data.mntd[, -1])

# Set them as a data frame, put back latitude and include categorical variables
data.mntd <- as.data.frame(data.mntd)
data.mntd$lat <- abs(data$lat)
data.mntd$lon <- data$lon
data.mntd$soil.fertility <- data$soil.fertility
data.mntd$soil.salinity <- data$soil.salinity

#======================================================================================

# Statistical analysis of phylogenetic alpha diversity

### SESmpd ###
# Empty model to start stepAIC for SESmpd
empty.model <- lm(sesmpd.qian ~ 1, data = data.mpd)

# Full model to finish stepAIC for SESmpd
full.model <- formula(lm(sesmpd.qian ~., data.mpd[, -7]))# Excluding longitude

# Stepwise selection based on AIC values
stepAIC(empty.model, direction = "both", scope = full.model)

# Final model SESmpd
final.model.mpd <- gls(sesmpd.qian ~ precip.driest.month + scale(lat) + soil.salinity, 
                       data = data.mpd,
                       correlation = corExp(form = ~lon + lat))


# Summary of the final model for SESmpd
summary(final.model.mpd)

# VIF of the final model for SESmpd
vif(final.model.mpd)

# Pseudo R square
cor(data.mpd$sesmpd.qian, predict(final.model.mpd))^2

# Plot the final model with standardised coefficients
summary.ses.mpd.plot <- plot_model(final.model.mpd, show.values  = TRUE, 
                                   group.terms = c(1,2,3), value.size = 5, wrap.labels = 22, vline.color = "light grey",
                                   title = "", colors = c("#009E73", "#009E73", "#666666"), axis.labels = c("Soil slightly saline", "Latitude", "Precipitation driest month"),
                                   axis.title = "Standardised coefficient")+
  theme_classic()+
  theme(axis.title.y = element_text(size = 20))+
  theme(axis.title.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  labs(title = expression("pseudo-R"^2 ~"= 0.31"))+
  theme(title = element_text(size = 20))


### SESmntd ###
# Empty model to start stepAIC for SESmntd
empty.model <- lm(sesmntd.qian ~ 1, data = data.mntd)

# Full model to finish stepAIC for SESmntd
full.model <- formula(lm(sesmntd.qian ~., data.mntd[, -7]))

# Stepwise selection based on AIC values
stepAIC(empty.model, direction = "both", scope = full.model)


# Final model SESmntd
final.model.mntd <- gls(sesmntd.qian ~ clim.stab.temp + soil.fertility + 
                          max.temp.warmest.month + soil.salinity + precip.driest.month, 
                        data = data.mntd,
                        correlation = corExp(form = ~lon + lat))

# Summary of the final model for SESmntd
summary(final.model.mntd)


#Pseudo R square
cor(data.mntd$sesmntd.qian,predict(final.model.mntd))^2


#vif of the final model for ses mntd
vif(final.model.mntd)


# Plot the final model with standardised coefficients
summary.ses.mntd.plot <- plot_model(final.model.mntd, show.values = T, value.offset = 0.3, vline = "light grey",
                                    title = "", order.terms = c(4,1,6,5,3,2), group.terms = c(4,1,6,5,3,2), value.size = 5, wrap.labels = 21,
                                    axis.labels = c("Dystrophic soil", "Mesotrophic soil", "Soil slightly saline", "Precipitation driest month", "Historical temperature instability", 
                                                    "Max. temp. warmest month"),
                                    colors = c("#666666", "#009E73", "#009E73", "#009E73", "#009E73", "#009E73"),
                                    axis.title = "Standardised coefficient")+
  theme_classic()+
  theme(axis.title.y = element_text(size = 20))+
  theme(axis.title.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.text.x = element_text(size=20))+
  labs(title = expression("pseudo-R"^2 ~"= 0.28"))+
  theme(title = element_text(size = 20))

#======================================================================================

#Plot main figure 2

# Prepare for high resolution figure
tiff("summary.gls.figures.new.tiff", units="in", width=15, height=7, res=600)

# Plot both figures together
ggarrange(summary.ses.mpd.plot, summary.ses.mntd.plot,
          labels = c("a)", "b)"), font.label = list(size = 20))

# Clean the plot environment
dev.off()



#======================================================================================

# Statistical analysis of phylogenetic beta diversity

# Get geographical coordinates
longlat <- cbind(data$lon, data$lat)

# Create geographic matrix
coord <- spDists(as.matrix(longlat), longlat = T, diagonal = F)
rownames(coord) <- data$code
colnames(coord) <- data$code

# Create geographic distance between assemblages
geo.dist <- as.dist(coord)


#===========================================================================================
# Principal component analysis

# PCA environment based on 15 environmental variables (columns 20 to 34 in the dataset)
pca.env <- prcomp(data[, 20:34], scale = T)
pca.env <- scores(pca.env, display = "sites", choice = 1:3)

# Create environmental distance between assemblages
env.dist <- dist(pca.env, method = "euclidian")

#===========================================================================================

# Create temperature distance between assemblages
temp.dist <- dist(data$annual.mean.temp, method = "euclidian")

#===========================================================================================

# Create precipitation distance between assemblages
precip.dist <- dist(data$annual.precip, method = "euclidian")

#===========================================================================================

# Extract residuals from a linear model with measures of phylobeta as 
# response variable and either geographic or environmental distance as predictors
# using both linear and quadratic terms.

# Residuals of phylobeta at deep phylogenetic levels
resid.dpw.space <- lm(phylobeta.dpw ~ poly(geo.dist, 2))$residuals
resid.dpw.env <- lm(phylobeta.dpw ~ poly(env.dist, 2))$residuals

# Residuals of phylobeta at shallow phylogenetic levels
resid.dnn.space <- lm(phylobeta.dnn ~ poly(geo.dist, 2))$residuals
resid.dnn.env <- lm(phylobeta.dnn ~ poly(env.dist, 2))$residuals

#Get residuals from annual mean temperature
resid.dpw.temp <- lm(phylobeta.dpw ~ poly(temp.dist, 2))$residuals
resid.dnn.temp <- lm(phylobeta.dnn ~ poly(temp.dist, 2))$residuals

#Get residuals from annual precipitation
resid.dpw.precip <- lm(phylobeta.dpw ~ poly(precip.dist, 2))$residuals
resid.dnn.precip <- lm(phylobeta.dnn ~ poly(precip.dist, 2))$residuals

#mantel test

#Mantel test using residuals of phylogenetic beta diversity, geographic and environmental distances.
ecodist::mantel(resid.dpw.space~ env.dist)
ecodist::mantel(resid.dpw.env~ geo.dist)
ecodist::mantel(resid.dnn.space~ env.dist)
ecodist::mantel(resid.dnn.env~ geo.dist)

#Mantel test using residuals of phylogenetic beta diversity, geographical distance and temperature distance.
ecodist::mantel(resid.dpw.space~ temp.dist)
ecodist::mantel(resid.dpw.temp~ geo.dist)
ecodist::mantel(resid.dnn.space~ temp.dist)
ecodist::mantel(resid.dnn.temp~ geo.dist)

#Mantel test using residuals of phylogenetic beta diversity, geographical distance and precipitation distance.
ecodist::mantel(resid.dpw.space~ precip.dist)
ecodist::mantel(resid.dpw.precip~ geo.dist)
ecodist::mantel(resid.dnn.space~ precip.dist)
ecodist::mantel(resid.dnn.precip~ geo.dist)

#===========================================================================================

# Get residuals from all predictors for phylobeta using Simpson pairwise dissimilarity index
resid.true.phylobeta.space <- lm(simpson.true.phylo.beta ~ poly(geo.dist, 2))$residuals
resid.true.phylobeta.env <- lm(simpson.true.phylo.beta ~ poly(env.dist, 2))$residuals
resid.true.phylobeta.temp <- lm(simpson.true.phylo.beta ~ poly(temp.dist, 2))$residuals
resid.true.phylobeta.precip <- lm(simpson.true.phylo.beta ~ poly(precip.dist, 2))$residuals

#Mantel tests
ecodist::mantel(resid.true.phylobeta.env ~ geo.dist)
ecodist::mantel(resid.true.phylobeta.temp ~ geo.dist)
ecodist::mantel(resid.true.phylobeta.precip ~ geo.dist)

ecodist::mantel(resid.true.phylobeta.space ~ env.dist)
ecodist::mantel(resid.true.phylobeta.space ~ temp.dist)
ecodist::mantel(resid.true.phylobeta.space ~ precip.dist)

#===========================================================================================

#Plot main figure 3

# Put names in two rows for the legend of the following figure
phylobeta.residuals <- gsub(" ", "\n", "Phylobeta (residuals)")


# Plot residuals of Dpw~geo.dist with environmental distance
plot.phylobeta.dpw.env <- ggplot(data = NULL, aes(env.dist, resid.dpw.space))+
  geom_point(aes(color = resid.dpw.space), alpha = 0.3)+
  geom_smooth(method = "loess", se = F, color = "dark blue", lwd = 1.2)+
  theme_classic()+
  scale_y_continuous(name="", breaks=pretty(resid.dpw.space, n = 6))+
  scale_x_continuous(name="")+
  theme(axis.title.x  = element_text(size = 20))+
  theme(axis.title.y  = element_text(size = 20))+
  scale_color_viridis_c(space = "Lab")+
  labs(color = phylobeta.residuals)

# Plot residuals of Dpw~env.dist with geographic distance
plot.phylobeta.dpw.space.env<- ggplot(data = NULL, aes(geo.dist, resid.dpw.env))+
  geom_point(aes(color = resid.dpw.env), alpha = 0.3)+
  geom_smooth(method = "loess", se = F, color = "dark blue", lwd = 1.2)+
  theme_classic()+
  scale_y_continuous(name="Dpw phylobeta diversity (residuals)", breaks=pretty(resid.dpw.env, n = 6))+
  scale_x_continuous(name="")+
  theme(axis.title.x  = element_text(size = 14))+
  theme(axis.title.y  = element_text(size = 14))+
  scale_color_viridis_c(space = "Lab")+
  labs(color = phylobeta.residuals)

# Plot residuals of Dnn~geo.dist with environmental distance
plot.phylobeta.dnn.env <- ggplot(data = NULL, aes(env.dist, resid.dnn.space))+
  geom_point(aes(color = resid.dnn.space), alpha = 0.3)+
  geom_smooth(method = "loess", se = F, color = "dark blue", lwd = 1.2)+
  theme_classic()+
  scale_y_continuous(name="", breaks=pretty(resid.dnn.space, n = 5))+
  scale_x_continuous(name="Environmental distance")+
  theme(axis.title.x  = element_text(size = 14))+
  theme(axis.title.y  = element_text(size = 14))+
  scale_color_viridis_c(space = "Lab")+
  labs(color = phylobeta.residuals)

# Plot residuals of Dnn~env.dist with geographic distance
plot.phylobeta.dnn.space.env<- ggplot(data = NULL, aes(geo.dist, resid.dnn.env))+
  geom_point(aes(color = resid.dnn.env), alpha = 0.3)+
  geom_smooth(method = "loess", se = F, color = "dark blue", lwd = 1.2)+
  theme_classic()+
  scale_y_continuous(name="Dnn phylobeta diversity (residuals)", breaks=pretty(resid.dnn.env, n = 6))+
  scale_x_continuous(name="Geographic distance (km)")+
  theme(axis.title.x  = element_text(size = 14))+
  theme(axis.title.y  = element_text(size = 14))+
  scale_color_viridis_c(space = "Lab")+
  labs(color = phylobeta.residuals)



# Prepare for high resolution figure
tiff("phylobeta.residuals.env.tiff", units="in", width=12, height=12, res=600)

# Plot all figures together
ggarrange(plot.phylobeta.dpw.space.env, plot.phylobeta.dpw.env, 
          plot.phylobeta.dnn.space.env, plot.phylobeta.dnn.env,
          labels = c("a)", "b)", "c)", "d)"),
          ncol = 2, nrow = 2)

# Close the plot window
dev.off()

#===========================================================================================

# Climatic instability plot (Appendix S1, figure S1.2)

# Prepare for high resolution figure
tiff("Climatic.instability.tiff", units="in", width=6, height=4, res=600)

# Plot the figure
ggplot(data, aes(abs(lat), data$clim.stab.temp))+
  geom_line(lwd = 1.5)+
  theme_classic()+
  scale_y_continuous(name = "Historical temperature instability (ºC)")+
  scale_x_continuous(name = "Latitude")+
  theme(axis.title.y = element_text(size = 15))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.text.x = element_text(size=12))

# Close the plot window
dev.off()
#===========================================================================================

# Precipitation plots (Appendix 1, figure S1.3)

# Precipitation seasonality
precip.seaso <- ggplot(data, aes(abs(lat), precip.seasonality))+
  geom_point(alpha = 0.3)+
  geom_smooth(color = "red", se = F)+
  theme_classic()+
  scale_y_continuous(name = "Precipitation seasonality (CV)")+
  scale_x_continuous(name = "Latitude (N or S)")+
  theme(axis.title.y = element_text(size = 15))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.text.x = element_text(size=12))

# Precipitation of the driest month
precip.driest <- ggplot(data, aes(abs(lat), precip.driest.month))+
  geom_point(alpha = 0.3)+
  geom_smooth(color = "red", se = F)+
  theme_classic()+
  scale_y_continuous(name = "Precipitation of the driest month")+
  scale_x_continuous(name = "Latitude (N or S)")+
  theme(axis.title.y = element_text(size = 15))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.text.x = element_text(size=12))

# High resolution figure
tiff("Precipitation.plots.tiff", units="in", width=9, height=4, res=600)

# Plot the figure
ggarrange(precip.seaso, precip.driest,
          labels = c("a)", "b)",
                     ncol = 2, nrow = 1))

# Close the plot window
dev.off()
#===========================================================================================

# Plot figure S2.1

# High resolution figure
tiff("most.abundant.families.tiff", units="in", width=14, height=7, res=600)

# Plot the figure. Data for this figure come from community data with
#one extra column (family).
ggplot(families, aes(abs(Lat), value))+
  geom_point()+
  geom_smooth(se = F, color = "red")+
  facet_wrap(~variable, nrow = 2, ncol = 5)+
  theme_bw()+
  scale_y_continuous(name = "Number of species / family")+
  scale_x_continuous(name = "Latitude")+
  theme(axis.title.y = element_text(size = 17))+
  theme(axis.title.x = element_text(size = 17))+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.text.x = element_text(size=12))+
  theme(strip.text = element_text(size = 15))

# Close the plot window
dev.off()

#===========================================================================================
# Plot figure S4.1

# Plot residuals of Dpw~geo.dist with temperature 
plot.phylobeta.dpw.temp <- ggplot(data = NULL, aes(temp.dist, resid.dpw.space))+
  geom_jitter(alpha = 0.3, width = 0.3)+
  geom_smooth(method = "loess", se = F, color = "red", lwd = 1.2)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_y_continuous(name="", breaks=pretty(resid.dpw.space, n = 6))+
  scale_x_continuous(name="")+
  theme(axis.title.x  = element_text(size = 14))+
  theme(axis.title.y  = element_text(size = 14))

# Plot residuals of Dpw~temp.dist with geographic distance
plot.phylobeta.dpw.space.temp<- ggplot(data = NULL, aes(geo.dist, resid.dpw.temp))+
  geom_point(alpha = 0.3)+
  geom_smooth(method = "loess", se = F, color = "red", lwd = 1.2)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_y_continuous(name="Dpw phylobeta diversity (residuals)", breaks=pretty(resid.dpw.temp, n = 6))+
  scale_x_continuous(name="")+
  theme(axis.title.x  = element_text(size = 14))+
  theme(axis.title.y  = element_text(size = 14))

# Plot residuals of Dnn~geo.dist with temperature distance
plot.phylobeta.dnn.temp <- ggplot(data = NULL, aes(temp.dist, resid.dnn.space))+
  geom_jitter(alpha = 0.3, width = 0.3)+
  geom_smooth(method = "loess", se = F, color = "red", lwd = 1.2)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_y_continuous(name="", breaks=pretty(resid.dnn.space, n = 6))+
  scale_x_continuous(name="Temperature distance (ºC)")+
  theme(axis.title.x  = element_text(size = 14))+
  theme(axis.title.y  = element_text(size = 14))

# Plot residuals of Dnn~temp.dist with geographic distance
plot.phylobeta.dnn.space.temp<- ggplot(data = NULL, aes(geo.dist, resid.dnn.temp))+
  geom_point(alpha = 0.3)+
  geom_smooth(method = "loess", se = F, color = "red", lwd = 1.2)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_y_continuous(name="Dnn phylobeta diversity (residuals)", breaks=pretty(resid.dnn.temp, n = 6))+
  scale_x_continuous(name="Geographic distance (km)")+
  theme(axis.title.x  = element_text(size = 14))+
  theme(axis.title.y  = element_text(size = 14))

# Prepare for high resolution figure
tiff("phylobeta.residuals.temperature.tiff", units="in", width=12, height=12, res=600)

ggarrange(plot.phylobeta.dpw.space.temp, plot.phylobeta.dpw.temp, 
          plot.phylobeta.dnn.space.temp, plot.phylobeta.dnn.temp,
          labels = c("a)", "b)", "c)", "d)"),
          ncol = 2, nrow = 2)

dev.off()
#===========================================================================================

# Plot figure S4.2

# Plot residuals of Dpw~geo.dist with precipitation distance
plot.phylobeta.dpw.precip <- ggplot(data = NULL, aes(precip.dist, resid.dpw.space))+
  geom_jitter(alpha = 0.3, width = 0.3)+
  geom_smooth(method = "loess", se = F, color = "red", lwd = 1.2)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_y_continuous(name="", breaks=pretty(resid.dpw.space, n = 6))+
  scale_x_continuous(name="")+
  theme(axis.title.x  = element_text(size = 14))+
  theme(axis.title.y  = element_text(size = 14))

# Plot residuals of Dpw~precip.dist with geographic distance
plot.phylobeta.dpw.space.precip<- ggplot(data = NULL, aes(geo.dist, resid.dpw.precip))+
  geom_point(alpha = 0.3)+
  geom_smooth(method = "loess", se = F, color = "red", lwd = 1.2)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_y_continuous(name="Dpw phylobeta diversity (residuals)", breaks=pretty(resid.dpw.precip, n = 6))+
  scale_x_continuous(name="")+
  theme(axis.title.x  = element_text(size = 14))+
  theme(axis.title.y  = element_text(size = 14))

# Plot residuals of Dnn~geo.dist with precipitation distance
plot.phylobeta.dnn.precip <- ggplot(data = NULL, aes(precip.dist, resid.dnn.space))+
  geom_jitter(alpha = 0.3, width = 0.3)+
  geom_smooth(method = "loess", se = F, color = "red", lwd = 1.2)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_y_continuous(name="", breaks=pretty(resid.dnn.space, n = 6))+
  scale_x_continuous(name="Precipitation distance (mm)")+
  theme(axis.title.x  = element_text(size = 14))+
  theme(axis.title.y  = element_text(size = 14))

# Plot residuals of Dnn~precip.dist with geographic distance
plot.phylobeta.dnn.space.precip<- ggplot(data = NULL, aes(geo.dist, resid.dnn.precip))+
  geom_point(alpha = 0.3)+
  geom_smooth(method = "loess", se = F, color = "red", lwd = 1.2)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_y_continuous(name="Dnn phylobeta diversity (residuals)", breaks=pretty(resid.dnn.precip, n = 6))+
  scale_x_continuous(name="Geographic distance (km)")+
  theme(axis.title.x  = element_text(size = 14))+
  theme(axis.title.y  = element_text(size = 14))

# Prepare for high resolution figure
tiff("phylobeta.residuals.precipitation.tiff", units="in", width=12, height=12, res=600)

# Plot all figures together
ggarrange(plot.phylobeta.dpw.space.precip, plot.phylobeta.dpw.precip, 
          plot.phylobeta.dnn.space.precip, plot.phylobeta.dnn.precip,
          labels = c("a)", "b)", "c)", "d)"),
          ncol = 2, nrow = 2)

# Close the plot window
dev.off()
#===========================================================================================

# Plot figure S4.3

# Plot residuals of phylo.simpson ~ env.dist with geographic distance
phylosimp.space.env <- ggplot(data = NULL, aes(geo.dist, resid.true.phylobeta.env))+
  geom_point(alpha = 0.3)+
  geom_smooth(method = "loess", se = F, color = "red", lwd = 1.2)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_y_continuous(name="Phylobeta Simpson (residuals)", breaks=pretty(resid.true.phylobeta.env, n = 6))+
  scale_x_continuous(name="Geographic distance (km)")+
  theme(axis.title.x  = element_text(size = 18))+
  theme(axis.title.y  = element_text(size = 18))+
  theme(axis.text=element_text(size=18))

# Plot residuals of phylo.simpson ~ temp.dist with geographic distance
phylosimp.space.temp <- ggplot(data = NULL, aes(geo.dist, resid.true.phylobeta.temp))+
  geom_point(alpha = 0.3)+
  geom_smooth(method = "loess", se = F, color = "red", lwd = 1.2)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_y_continuous(name="", breaks=pretty(resid.true.phylobeta.temp, n = 6))+
  scale_x_continuous(name="Geographic distance (km)")+
  theme(axis.title.x  = element_text(size = 18))+
  theme(axis.title.y  = element_text(size = 18))+
  theme(axis.text=element_text(size=18))

# Plot residuals of phylo.simpson ~ precip.dist with geographic distance
phylosimp.space.precip <- ggplot(data = NULL, aes(geo.dist, resid.true.phylobeta.precip))+
  geom_point(alpha = 0.3)+
  geom_smooth(method = "loess", se = F, color = "red", lwd = 1.2)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_y_continuous(name="", breaks=pretty(resid.true.phylobeta.precip, n = 6))+
  scale_x_continuous(name="Geographic distance (km)")+
  theme(axis.title.x  = element_text(size = 18))+
  theme(axis.title.y  = element_text(size = 18))+
  theme(axis.text=element_text(size=18))


# Plot residuals of phylo.simpson ~ geo.dist with environmental distance
phylosimp.env.space <- ggplot(data = NULL, aes(env.dist, resid.true.phylobeta.space))+
  geom_point(alpha = 0.3)+
  geom_smooth(method = "loess", se = F, color = "red", lwd = 1.2)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_y_continuous(name="Phylobeta Simpson (residuals)", breaks=pretty(resid.true.phylobeta.space, n = 6))+
  scale_x_continuous(name="Environmental distance")+
  theme(axis.title.x  = element_text(size = 18))+
  theme(axis.title.y  = element_text(size = 18))+
  theme(axis.text=element_text(size=18))

# Plot residuals of phylo.simpson ~ geo.dist with temperature distance
phylosimp.temp.space <- ggplot(data = NULL, aes(temp.dist, resid.true.phylobeta.space))+
  geom_point(alpha = 0.3)+
  geom_smooth(method = "loess", se = F, color = "red", lwd = 1.2)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_y_continuous(name="", breaks=pretty(resid.true.phylobeta.space, n = 6))+
  scale_x_continuous(name="Temperature distance (ºC)")+
  theme(axis.title.x  = element_text(size = 18))+
  theme(axis.title.y  = element_text(size = 18))+
  theme(axis.text=element_text(size=18))

# Plot residuals of phylo.simpson ~ geo.dist with precipitation distance
phylosimp.precip.space <- ggplot(data = NULL, aes(precip.dist, resid.true.phylobeta.space))+
  geom_point(alpha = 0.3)+
  geom_smooth(method = "loess", se = F, color = "red", lwd = 1.2)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_y_continuous(name="", breaks=pretty(resid.true.phylobeta.space, n = 6))+
  scale_x_continuous(name="Precipitation distance (mm)")+
  theme(axis.title.x  = element_text(size = 18))+
  theme(axis.title.y  = element_text(size = 18))+
  theme(axis.text=element_text(size=18))


# Prepare for high resolution figure
tiff("phylobeta.simpson.tiff", units="in", width=18, height=12, res=600)

# Plot the figure
ggarrange(phylosimp.space.env, phylosimp.space.temp, phylosimp.space.precip,
          phylosimp.env.space, phylosimp.temp.space,phylosimp.precip.space,
          labels = c("a)", "b)", "c)", "d)", "e)", "f)"),
          ncol = 3, nrow = 2, font.label = list(size = 18))

# Close the plot window
dev.off()

# END

