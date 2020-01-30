
# Part 3, which includes the figures in the main text.

# Load required packages

library(openxlsx)
library(here)
library(tidyverse)
library(mapdata)
library(spdep)
library(geosphere)
library(raster)
library(gridExtra)
library(grid)
library(ggsn)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggpubr)


#======================================================================================

# Brazilian map with restinga assemblages (main figure 1)

# The map included in this study was created following a tutorial available at
# https://rstudio-pubs-static.s3.amazonaws.com/176768_ec7fb4801e3a4772886d61e65885fbdd.html
# You can find the layer "Brazil_biomes" from there.

# Read shapefile (Brazil_biomes layer)
br <- rgdal::readOGR(dsn = "C:\\Users\\jcmas\\Nextcloud\\PhD_UNITARTU\\Current projects\\Phylogenetic-diversity-Restingas\\Brazil_biomes", layer = "Brazil_biomes")

# Load the full dataset to get geographic coordinates and ecoregions
data <- openxlsx::read.xlsx(here("Data", "data.full.xlsx"))

#======================================================================================

# Check projection
proj4string(br)

# Check names
names(br)
head(br$name)

# Select atlantic biome                     
atlantic <- br[br@data$name=="Mata AtlÃ¢ntica",]
MA <- atlantic
MAf <- fortify(MA)#Fortify to bw used in ggplot2

# Plot atlantic biome in ggplot2
brazil.map <- borders("worldHires", regions = "Brazil", fill = "white", colour = "black")
bra.map <- map_data("worldHires","Brazil")

#The map has been limited to include longitude between -60 and -34, and latitude between 0 and -35.
map.1 <- ggplot()+
  brazil.map+
  coord_equal()+
  coord_map(xlim=c(-60, -34), ylim=c(-35, 0))


# Rename ecoregion column from the dataset and transform it into factor
data$eco.region <- as.factor(data$eco.region)

levels(data$eco.region)[match("Northeastern Atlantic Shorelands",levels(data$eco.region))] <- "Northern restingas"
levels(data$eco.region)[match("Central Atlantic Shorelands",levels(data$eco.region))] <- "Central restingas"
levels(data$eco.region)[match("Southern Atlantic Shorelands",levels(data$eco.region))] <- "Southern restingas"
data$eco.region <- factor(data$eco.region, levels = c("Northern restingas", "Central restingas", "Southern restingas"))

# Fortify data.frame with geographic coordinates of sites (taken from the dataset)
coord <- data[, 8:9]
coord <- fortify(coord)

# Include ecoregions in the coord data frame to be plotted.
coord$eco.region <- factor(data$eco.region, levels = c("Northern restingas", 
                                                       "Central restingas",
                                                       "Southern restingas"))
# Plot brazil map with the Atlantic forest domain 
map.2 <- map.1 +
  geom_polygon(data = MAf, aes(x = long, y = lat, group = group), fill = "grey", alpha = 0.5)+
  theme_bw()

#======================================================================================

# Get South American map

# Global map
world <- map_data("world")

# It needs a vector with all countries of South America as an input
countrylist <- read.table(here("Data", "countrylist.txt"))

# Get coordinates of countries
lat <- map_data("world",region=countrylist$V1)

# Remove some islands for visual effects
lat <- lat %>% 
  filter(long >= -82)

idx <- which(lat$long<=-75 & lat$lat <= -20)

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
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#======================================================================================

# Plot the Brazilian map with the Atlantic forest domain and assemblages used in this study

# Plot the map
g2 <- map.2+
  geom_point(data = coord, aes(lon, lat, color = eco.region), size = 1.15, alpha = 0.8)+
  theme_bw()+
  labs(color = "")+
  scale_colour_viridis_d()+
  ggsn::scalebar(data = NULL, dist = 500, dist_unit = "km", model = 'WGS84',
                 transform = T, st.size=2.5, height=0.09, st.dist = 0.05,
                 x.min = -45,
                 x.max = -37,
                 y.min = -34,
                 y.max = -25)+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(legend.text = element_text(size = 17))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_x_continuous(name = 'Longitude', labels = c("60ºW", "50ºW", "40ºW"), breaks = c(-60, -50, -40))+
  scale_y_continuous(name = 'Latitude', labels =  c("0º", "10ºS", "20ºS", "30ºS"), breaks = c(0, -10, -20, -30))+
  theme(axis.title.y = element_text(size = 17))+
  theme(axis.title.x = element_text(size = 17))+
  geom_segment(aes(x = -36, y = -4, xend = -36, yend = -2), arrow = arrow(length = unit(0.02, "npc")))+
  annotate("text", x = -36, y = -1, label = "N")

# Figure 1

# Plot both maps together (South American map as an inset map) and export them as a tiff file
tiff(file="map.tiff", units = "in", w= 7, h= 6, res=600)#It could be saved with ggsave as well

grid.newpage()
main <- viewport(width = 1, height = 1.3, x = 0.5, y = 0.5)# Plot area for the main map
print(g2, vp = Main)# Main map 

inset <- viewport(width = 0.28, height = 0.3, x = 0.21, y = 0.75)# Plot area for the inset map
print(g1, vp = inset)# South American map

# Close the plot window
dev.off()

#======================================================================================

# Figure 2

# Plot the relationship between species richness and latitude
# NOTE: It need the object "pred" from part 2 (statistical analysis)
# to plot the 95% confidence intervals from the model.

richness.plot <- ggplot(data, aes(abs(lat), species))+
  geom_point(aes(color = data$eco.region))+
  geom_line(data=pred, color = "dark blue", lwd = 1.2, alpha = 0.5)+
  geom_ribbon(data=pred,aes(ymin=lwr,ymax=upr), alpha=0.2)+
  scale_x_continuous(name = "Latitude")+
  scale_y_continuous(name = "Species richness")+
  theme(axis.title.y = element_text(size = 17))+
  theme(axis.title.x = element_text(size = 17))+
  scale_color_viridis_d()+
  theme_bw()+
  labs(color = "")+
  guides(colour = guide_legend(override.aes = list(size=4)))

#Save it
ggsave("richness.plot.tiff", richness.plot, width = 5, height = 3)

#======================================================================================

# Figure 3

# Because the 95% confidence intervals for each relationship
# were calculated separately from each model, I could not figure out
# how to plot them using "face_grid" or "facet_wrap" from ggplot2, which would be much easier
# Perhaps using some gather or spread functions, but I did not try it. Therefore, each figure
# was produced separately using the output of each model to plot lines and their 95%CI.
# and they are stored in a different folder called Analysis_family.

#======================================================================================

# Figure 4

# Plot the standardised coefficients and 95%CI related to both indexes
# of phylogenetic alpha diversity

# SESmpd
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

# SESmntd
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


# Prepare for high resolution figure
tiff("summary.gls.figures.new.tiff", units="in", width=15, height=7, res=600)

# Plot both figures together
ggarrange(summary.ses.mpd.plot, summary.ses.mntd.plot,
          labels = c("a)", "b)"), font.label = list(size = 10))

# Clean the plot environment
dev.off()

#======================================================================================

# Figure 5

# It requires the output of the linear regression between phylobeta diversity
# and predictors (resid.dpw.space, resid.dpw.env, resid.dnn.space, resid.dnn.env), 
# and the geographic or environmental dissimilarity distances (all calculated the part 2)


# Put names in two rows for the legend of the following figure
phylobeta.residuals <- gsub(" ", "\n", "Phylobeta (residuals)")


# Plot residuals of Dpw~geo.dist with environmental distance
plot.phylobeta.dpw.env <- ggplot(data = NULL, aes(env.dist, resid.dpw.space))+
  geom_point(aes(color = resid.dpw.space), alpha = 0.3)+
  theme_bw()+
  scale_y_continuous(name="", breaks=pretty(resid.dpw.space, n = 6))+
  scale_x_continuous(name="")+
  theme(axis.title.x  = element_text(size = 20))+
  theme(axis.title.y  = element_text(size = 20))+
  scale_color_viridis_c(space = "Lab")+
  labs(color = phylobeta.residuals)+
  annotate("text", x = 0.5, y = 5.2, label = "n.s.")

# Plot residuals of Dpw~env.dist with geographic distance
plot.phylobeta.dpw.space.env <- ggplot(data = NULL, aes(geo.dist, resid.dpw.env))+
  geom_point(aes(color = resid.dpw.env), alpha = 0.3)+
  theme_bw()+
  scale_y_continuous(name="Dpw (residuals)", breaks=pretty(resid.dpw.env, n = 6))+
  scale_x_continuous(name="")+
  theme(axis.title.x  = element_text(size = 14))+
  theme(axis.title.y  = element_text(size = 14))+
  scale_color_viridis_c(space = "Lab")+
  labs(color = phylobeta.residuals)+
  annotate("text", x = 100, y = 5.2, label = "n.s.")

# Plot residuals of Dnn~geo.dist with environmental distance
plot.phylobeta.dnn.env <- ggplot(data = NULL, aes(env.dist, resid.dnn.space))+
  geom_point(aes(color = resid.dnn.space), alpha = 0.3)+
  geom_line(stat="smooth", method = "lm", se = F, color = "blue", lwd = 1.1, alpha = 0.5)+
  theme_bw()+
  scale_y_continuous(name="", breaks=pretty(resid.dnn.space, n = 5))+
  scale_x_continuous(name="Environmental distance")+
  theme(axis.title.x  = element_text(size = 14))+
  theme(axis.title.y  = element_text(size = 14))+
  scale_color_viridis_c(space = "Lab")+
  labs(color = phylobeta.residuals)+
  annotate("text", x = 12, y = -15, label = "r = 0.32 (p = 0.001)")


# Plot residuals of Dnn~env.dist with geographic distance
plot.phylobeta.dnn.space.env <- ggplot(data = NULL, aes(geo.dist, resid.dnn.env))+
  geom_point(aes(color = resid.dnn.env), alpha = 0.3)+
  geom_line(stat="smooth", method = "lm", se = F, color = "blue", lwd = 1.1, alpha = 0.5)+
  theme_bw()+
  scale_y_continuous(name="Dnn (residuals)", breaks=pretty(resid.dnn.env, n = 6))+
  scale_x_continuous(name="Geographic distance (km)")+
  theme(axis.title.x  = element_text(size = 14))+
  theme(axis.title.y  = element_text(size = 14))+
  scale_color_viridis_c(space = "Lab")+
  labs(color = phylobeta.residuals)+
  annotate("text", x = 2500, y = -30, label = "r = 0.25 (p = 0.001)")

# Prepare for high resolution figure
tiff("phylobeta.residuals.env.tiff", units="in", width=10, height=8, res=600)

# Plot all figures together
ggarrange(plot.phylobeta.dpw.space.env, plot.phylobeta.dpw.env, 
          plot.phylobeta.dnn.space.env, plot.phylobeta.dnn.env,
          labels = c("a)", "b)", "c)", "d)"),
          ncol = 2, nrow = 2)

# Close the plot window
dev.off()

# END
