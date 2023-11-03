library(ncdf4)
library(ggplot2)

# Import the data
ncpath <- "~/"
ncname <- "Bayelva_snow_geo"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")

# Open data 
ncin <- nc_open(ncfname)




# Create a data frame with the  first year of snow depth (ds)
df <- data.frame(LONG = ncvar_get(ncin,"lons"),
                 LAT = ncvar_get(ncin,"lats"),
                 #elev = ncvar_get(ncin,"zs"),
                 DS = ncvar_get(ncin,"ds")[,1],
                 TPI = ncvar_get(ncin,"TPI5s"),
                 DS_log = log(ncvar_get(ncin,"ds")[,1]))

## Ready for modelling!
TPI = ncvar_get(ncin,"TPI5s")

# look at the range of the data
initial_range <- diff(range(df[, "LONG"])) / 5
max_edge <- initial_range / 4
offset <- initial_range/2

mesh <- fmesher::fm_mesh_2d(
  loc = location_data[, c("LONG", "LAT")],
  max.edge = c(1, 2) * max_edge,
  offset = c(1,1.5)*offset,
  cutoff = max_edge/10 
)

fdmr::plot_mesh(mesh,spatial_data = location_data)

## Start defining the priors
# the prior range is the distance that the process should stop effecting, so in this case it is currently 20km away from the node center

spde1 <- INLA::inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(initial_range, 0.5),
  prior.sigma = c(1, 0.01)
)

spde2 <- INLA::inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(initial_range*5, 0.5),
  prior.sigma = c(1, 0.01)
)

## Creating our data frame into a SpatialPointsDataFrame ()

sp::coordinates(df) <- c("LONG", "LAT")

formula1 <- DS ~ 0 + Intercept(1)+
  f(
    main = df,
    model = spde1
  )


formula2 <- DS ~ 0 + Intercept(1)   +
  f(
    main = df,
    model = spde2
  )


inlabru_model1 <- inlabru::bru(formula1,
                               data = df,
                               family = "gaussian",
                               options = list(
                                 verbose = FALSE
                               )
)

inlabru_model2 <- inlabru::bru(formula2,
                               data = df,
                               family = "gaussian",
                               options = list(
                                 verbose = FALSE
                               )
)

summary(inlabru_model1)
summary(inlabru_model2)


plot(exp(inlabru_model1$summary.fitted.values$mean[1:nrow(df@data)]),exp(df@data$DS))
plot(exp(inlabru_model2$summary.fitted.values$mean[1:nrow(df@data)]),exp(df@data$DS))


plot(exp(inlabru_model1$summary.fitted.values$mean[1:nrow(df@data)]),exp(df@data$DS))
plot(exp(inlabru_model2$summary.fitted.values$mean[1:nrow(df@data)]),exp(df@data$DS))


fdmr::model_viewer(
  model_output = inlabru_model1,
  mesh = mesh,
  measurement_data = df,
  data_distribution = "Gaussian"
)

fdmr::model_viewer(
  model_output = inlabru_model2,
  mesh = mesh,
  measurement_data = df,
  data_distribution = "Gaussian"
)

