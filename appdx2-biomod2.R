rm(list=ls())
library(biomod2)
library(ggplot2)
library(gridExtra)
library(raster)
library(rasterVis)
library(mapview)
library(sf)

##########load in and clean data
occ <- read.csv("generalist_rare.csv")
names(occ)[names(occ) == '?..generalist'] <- 'generalist'
mapview(occ, xcol = "POINT_X", ycol = "POINT_Y", crs = 4269, grid = FALSE)

env <-
  raster::stack(c(
    ndvi = 'asc_ndvi_proj.asc',
    landuse = 'asc_landuse_proj.asc',
    slope = 'asc_slope_proj.asc',
    elev = 'asc_elev_proj.asc',
    hedge = 'asc_hedge.asc'
  ))

MyRespName <- 'generalist'
MyResp <- as.numeric(occ[,MyRespName])
myRespXY <- occ[,c("POINT_X", "POINT_Y")]

data <- 
  BIOMOD_FormatingData(
    resp.var = MyResp,  
    expl.var = env, 
    resp.xy =myRespXY,
    resp.name = "generalist",
    PA.nb.rep = 2,
    PA.nb.absences = 10000,
    PA.strategy = 'random'
  )
data
##########set modelling options
gen_options <-
  BIOMOD_ModelingOptions(
    GLM = list(type = 'quadratic', interaction.level = 1),
    GBM = list(n.trees = 1000),
    GAM = list(algo = 'GAM_mgcv'),
    MAXENT.Phillips = list(path_to_maxent.jar = "C:/Users/orie4396/Desktop/biomod",
                           linear = TRUE, 
                           quadratic = FALSE, 
                           product = FALSE, 
                           threshold = FALSE, 
                           hinge = FALSE))

gen_models <-
  BIOMOD_Modeling(
    data = data,
    models = c("MAXENT.Phillips","GBM", "RF", "GAM", "GLM"),
    models.options = ProLau_opt,
    NbRunEval = 10,
    DataSplit = 80,
    VarImport = 3,
    do.full.models = TRUE,
    modeling.id = "demo1"
  )

gen_models_scores <- get_evaluations(gen_models, as.data.frame=TRUE)

View(gen_models_scores)
write.csv(gen_models_scores, "generalist_evalscores.csv")

##########plot models evaluation scores
models_scores_graph(
  gen_models,
  by = "models",
  metrics = c("ROC", "TSS"),
  xlim = c(0.5,1),
  ylim = c(0.5,1)
)
models_scores_graph(
  gen_models,
  by = "cv_run",
  metrics = c("ROC", "TSS"),
  xlim = c(0.5,1),
  ylim = c(0.5,1)
)
models_scores_graph(
  gen_models,
  by = "data_set",
  metrics = c("ROC", "TSS"),
  xlim = c(0.5,1),
  ylim = c(0.5,1)
)

##########check variable importance, values closest to 1 are more important
(gen_models_var_import <- get_variables_importance(gen_models))
#make the mean of variable importance by algorithm
apply(gen_models_var_import, c(1,2), mean)
##########individual models response plots
gen_glm <- BIOMOD_LoadModels(gen_models, models = "GLM")
gen_gbm <- BIOMOD_LoadModels(gen_models, models = "GBM")
gen_rf <- BIOMOD_LoadModels(gen_models, models = "RF")
gen_gam <- BIOMOD_LoadModels(gen_models, models = "GAM")
gen_maxent <- BIOMOD_LoadModels(gen_models, models = "MAXENT.Phillips")

glm_eval_strip <- 
  biomod2::response.plot2(
    models = gen_glm,
    Data = get_formal_data(gen_models, 'expl.var'),
    show.variables= get_formal_data(gen_models, 'expl.var.names'), 
    do.bivariate=FALSE,
    fixed.var.metric = 'median',
    legend = FALSE,
    display_title = FALSE,
    data_species = get_formal_data(gen_models, 'resp.var')
  )
gbm_eval_strip <- 
  biomod2::response.plot2(
    models = gen_gbm,
    Data = get_formal_data(gen_models, 'expl.var'),
    show.variables= get_formal_data(gen_models, 'expl.var.names'), 
    do.bivariate=FALSE,
    fixed.var.metric = 'median',
    legend = FALSE,
    display_title = FALSE,
    data_species = get_formal_data(gen_models, 'resp.var')
  )
rf_eval_strip <- 
  biomod2::response.plot2(
    models = gen_rf,
    Data = get_formal_data(ProLau_models, 'expl.var'),
    show.variables= get_formal_data(gen_models, 'expl.var.names'), 
    do.bivariate=FALSE,
    fixed.var.metric = 'median',
    legend = FALSE,
    display_title = FALSE,
    data_species = get_formal_data(gen_models, 'resp.var')
  )
gam_eval_strip <- 
  biomod2::response.plot2(
    models = gen_gam,
    Data = get_formal_data(ProLau_models, 'expl.var'),
    show.variables= get_formal_data(gen_models, 'expl.var.names'), 
    do.bivariate=FALSE,
    fixed.var.metric = 'median',
    legend = FALSE,
    display_title = FALSE,
    data_species = get_formal_data(gen_models, 'resp.var')
  )
maxent_eval_strip <- 
  biomod2::response.plot2(
    models = gen_maxent,
    Data = get_formal_data(gen_models, 'expl.var'),
    show.variables= get_formal_data(gen_models, 'expl.var.names'), 
    do.bivariate=FALSE,
    fixed.var.metric = 'median',
    legend = FALSE,
    display_title = FALSE,
    data_species = get_formal_data(gen_models, 'resp.var')
  )

##########run the ensemble models
gen_ensemble_models <- 
  BIOMOD_EnsembleModeling(
    modeling.output = gen_models,
    em.by = 'all',
    eval.metric = 'TSS', 
    eval.metric.quality.threshold = 0.6, 
    models.eval.meth = c('TSS', 'ROC'),
    prob.mean = FALSE,
    prob.cv = TRUE, 
    committee.averaging = TRUE, 
    prob.mean.weight = TRUE,
    VarImport = 0
  )
##########assess ensemble models quality
(gen_ensemble_models_scores <- get_evaluations(gen_ensemble_models))
ProLau_ensemble_models_scores
#do model projections
gen_models_proj_current <- 
  BIOMOD_Projection(
    modeling.output = gen_models,
    new.env = env,
    proj.name = "current",
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )
gen_ensemble_models_proj_current<-
  BIOMOD_EnsembleForecasting(
    EM.output = gen_ensemble_models,
    projection.output = gen_models_proj_current,
    binary.meth = "TSS", 
    output.format = "img",
    do.stack = FALSE
  )

ensemblescores <- get_evaluations(gen_ensemble_models_proj_current)

