# Part 2 (Statistical analysis)

#======================================================================================
# Load required packages

library(openxlsx)
library(tidyverse)
library(here)
library(corrplot)
library(nlme)
library(MuMIn)
library(MASS)
library(ecodist)
library(car)
library(sp)
library(vegan)

#======================================================================================

#Load dataset
data <- openxlsx::read.xlsx(here("Data", "data.full.xlsx"))
str(data)

#======================================================================================
# Species richness analysis

model.richness <- gls(species ~ poly(abs(lat), 2), 
                      correlation = corExp(form = ~lon + lat),
                      data = data)

#Get confidence intervals for the model
#from http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#standard-errors-of-variance-estimates
#and https://stackoverflow.com/questions/14033551/r-plotting-confidence-bands-with-ggplot


#get predicted values
fit <- predict(model.richness)

#variace-covariance matrix
V <- vcov(model.richness)

#model matrix corresponding to the new data
X <- model.matrix( ~poly(abs(lat), 2), data = data)

#Standard errors
se.fit <- sqrt(diag(X %*% V %*% t(X)))

#data frame for the confidence intervals
pred <- with(data, data.frame(lat,
                              species = fit, lwr = fit-1.96*se.fit, upr = fit+1.96*se.fit))


#======================================================================================
#Preparation of predictors for analysis of phylogenetic alpha diversity

# Rename soil.fertility and transform it into factor
data$soil.fertility <- factor(data$soil.fertility)
levels(data$soil.fertility)[match("Dystrophic (10-25% TBS)",levels(data$soil.fertility))] <- "Dystrophic"
levels(data$soil.fertility)[match("Hypo",levels(data$soil.fertility))] <- "Hypo-dystrophic"
levels(data$soil.fertility)[match("Mesotrophic (25-50% TBS)",levels(data$soil.fertility))] <- "Mesotrophic"
data$soil.fertility <- factor(data$soil.fertility, levels = c("Hypo-dystrophic", "Dystrophic", "Mesotrophic"), ordered = F)

# Rename soil.salinity and transform it into factor
data$soil.salinity <- factor(data$soil.salinity, levels = c("Low.salinity", "Slightly.salinity"), ordered = F)
levels(data$soil.salinity)[match("Low.salinity",levels(data$soil.salinity))] <- "Low salinity"
levels(data$soil.salinity)[match("Slightly.salinity",levels(data$soil.salinity))] <- "Slightly salinity"


# Variables to check the correlation
vars.mpd <- c("lat", "annual.mean.temp", "temp.seasonality", "annual.precip", "precip.seasonality", "max.temp.warmest.month",
              "precip.driest.month", "clim.stab.temp", "clim.stab.precip")

# Perform pairwise correlation
data.mpd <- data[, vars.mpd]
corr.mpd <- cor(data.mpd)

# Plot the figure
corrplot(corr.mpd, type = "upper", diag = F, method = "number",
         number.cex = 1.5, number.font = 2, tl.cex = 1.5, tl.col = "black")

# Close the plot window 
dev.off()

#Excluded variables strongly correlated with latitude and limiting factors(r > 0.60) . It included
#annual mean temperature, temp.seasonality, precip seasonality, annual precipitation. 
vars.mpd <- c("lat", "max.temp.warmest.month", 
              "precip.driest.month", "clim.stab.temp", "clim.stab.precip", "sesmpd.qian")


# Scale the variables, except latitude (will be scaled after get absolute its values)
data.mpd <- data[, vars.mpd]
data.mpd <- scale(data.mpd[, -1])

# Put latitude back and include categorical variables (soil fertility and soil salinity)
data.mpd <- as.data.frame(data.mpd)
data.mpd$lat <- abs(data$lat)
data.mpd$lon <- data$lon
data.mpd$soil.fertility <- data$soil.fertility
data.mpd$soil.salinity <- data$soil.salinity

#======================================================================================
# Repeat the procedure for SESmntd

vars.mntd <- c("lat", "max.temp.warmest.month", 
               "precip.driest.month", "clim.stab.temp", "clim.stab.precip", "sesmntd.qian")

# Get the variables from the data frame and scale them
# except latitude (will be scaled after get absolute its values)
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
                       correlation = corExp(form = ~ lon + lat))

# Summary of the final model for SESmpd
summary(final.model.mpd)

# VIF of the final model for SESmpd
vif(final.model.mpd)

# Pseudo R square
cor(data.mpd$sesmpd.qian, predict(final.model.mpd))^2

#======================================================================================


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

# Extract residuals from a linear model with measures of phylobeta as 
# response variable and either geographic or environmental distance as predictors
# using both linear and quadratic terms.


# Residuals of phylobeta at deep phylogenetic levels
resid.dpw.space <- lm(phylobeta.dpw ~ poly(geo.dist, 2))$residuals
resid.dpw.env <- lm(phylobeta.dpw ~ poly(env.dist, 2))$residuals

# Residuals of phylobeta at shallow phylogenetic levels
resid.dnn.space <- lm(phylobeta.dnn ~ poly(geo.dist, 2))$residuals
resid.dnn.env <- lm(phylobeta.dnn ~ poly(env.dist, 2))$residuals

#mantel test

#Mantel test using residuals of phylobeta diversity, geographic and environmental distances.
ecodist::mantel(resid.dpw.space~ env.dist)
ecodist::mantel(resid.dpw.env~ geo.dist)
ecodist::mantel(resid.dnn.space~ env.dist)
ecodist::mantel(resid.dnn.env~ geo.dist)

# END
