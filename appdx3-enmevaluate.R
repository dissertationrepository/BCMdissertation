rm(list=ls())
library(ENMeval)
library(raster)
library(MASS)
-----

slope <- raster("asc_slope_proj.asc")
water <- raster("asc_water_proj.asc")
wood <- raster("asc_wood_proj.asc")
elev <- raster("asc_elev_proj.asc")
ndvi <- raster("asc_ndvi_proj.asc")
landuse <- raster("asc_landuse_proj.asc")
hedge <- raster("asc_hedge.asc")
env <- stack(slope, elev, aspect, ndvi, landuse, hedge)

occ <- read.csv("wood_rare.csv")[,-1]
View(occ)

###########make a bias file
occur.ras <- rasterize(occ, env, 1)
plot(occur.ras)

presences <- which(values(occur.ras) == 1)
pres.locs <- coordinates(occur.ras)[presences, ]

dens <- kde2d(pres.locs[,1], pres.locs[,2], n = c(nrow(occur.ras), ncol(occur.ras)), lims = c(extent(env)[1], extent(env)[2], extent(env)[3], extent(env)[4]))
dens.ras <- raster(dens, env)
dens.ras2 <- resample(dens.ras, env)
plot(dens.ras2)

writeRaster(dens.ras2, "wood_biasfile.asc", overwrite = TRUE)

#check how many potential background points you have available and select background points
length(which(!is.na(values(subset(env, 1)))))

bg <- xyFromCell(dens.ras2, sample(which(!is.na(values(subset(env, 1)))), 10000, prob=values(dens.ras2)[!is.na(values(subset(env, 1)))]))
colnames(bg) <- colnames(occ)
###########finetune evaluation
enmeval_results <- ENMevaluate(occ, env, bg, tune.args = list(fc = c("L","LQ","H", "LQH", "LQHP", "LQHPT"), rm = 1:5), partitions = "randomkfold", partition.settings = list(kfolds = 10), algorithm = "maxnet")
write.csv(enmeval_results@results, "urban_enmeval_results.csv")
View(enmeval_results)
