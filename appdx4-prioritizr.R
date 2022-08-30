# Code based on prioritizr.net/articles/prioritizr.html
rm(list=ls())
install.packages("c:/gurobi952/win64/R/gurobi_9.5-2.zip", repos = NULL)
library(prioritizr)
library(prioritizrdata)
library(raster)
library(rgdal)
library(sf)
library(dplyr)
library(sp)
library(slam)
library(gurobi)
options(tibble.width = Inf)

setwd("D:/final/prioritizr")
###########Input planning unit file
###each cell including connectivity cost value and a binary value for whether core area exists
pu_unproj <- st_read("aq_pu200inv2.shp")
names(pu_unproj)[names(pu_unproj) == 'grid_code'] <- 'Cost'

#modifying in binary column (logical vector) for protected area (labelled as TRUE)
pu_unproj <- pu_unproj %>%
  mutate(locked_in = (Field ==1))
is.logical(pu_unproj$locked_in)
View(pu_unproj)

plot(st_as_sf(pu_unproj[, "Cost"]), main = "Planning unit costs")
plot(st_as_sf(pu_unproj[, "locked_in"]), main = "Protected area coverage", pal = c("grey90", "darkgreen"))

env <-
  raster::stack(c(
    wood = 'prjwood_200.asc',
    water = 'prjaqua_200.asc'
  ))
env

#project planning unit shapefile to the same CRS as environmental variables
rst <- raster("prjaqua_200.asc")
pu <- st_transform(pu_unproj, crs(rst))
compareCRS(env, pu)

###########setting problem and solution
targets <- c(0.50)
targets <- c(0.30)
p1 <- problem(pu, env, cost_column= "Cost") %>%
  add_min_set_objective() %>%
  add_boundary_penalties(penalty = 0.005) %>%
  add_relative_targets(targets) %>%
  add_locked_in_constraints("locked_in") %>%
  add_binary_decisions()

print(p1)

s1 <- solve(p1)
plot(st_as_sf(s1[, "solution_1"]), main = "Prioritization",
     pal = c("grey90", "darkgreen"))

print(attr(s1, "objective"))
print(attr(s1, "runtime"))
print(attr(s1, "status"))
print(eval_cost_summary(p1, s1[, "solution_1"]), width = Inf)
s1
###########calculate how well are features represented in the solution 
pu_unproj$pa <- round(pu_unproj$locked_in)

# calculate feature representation statistics based on existing protected areas
tc_pa <- eval_target_coverage_summary(p1, pu[, "Field"])
print(tc_pa)
#View(tc_pa)
# calculate  feature representation statistics based on the prioritization
tc_s1 <- eval_target_coverage_summary(p1, s1[, "solution_1"])
print(tc_s1)
View(tc_s1)

###########calculating Ferrier's score (irreplaceability)
irrep_s1 <- eval_ferrier_importance(p1, s1["solution_1"])
print(irrep_s1)
# manually coerce values for planning units not selected in prioritization
# to NA, so that they are shown in white
irrep_s1$plot_total <- irrep_s1$total
irrep_s1$plot_total[s1$solution_1 < 0.5] <- NA_real_

# plot map of overall importance scores
plot(st_as_sf(irrep_s1[, "plot_total"]), main = "Overall importance")
