# pack_in_netcdf
R script to convert non-regularly gridded files into netcdf 

I am going to use some libraries available in R to help on this task

````R

library(ncdf4)
library(raster)
library(akima)
library(RColorBrewer)
spec=brewer.pal(n=11, "RdYlBu")
````

We are loading our initial data which is a .Rdata formatted saved workspace from a previous step. 
They are triangularly gridded dataframes containing X and Y, time and out variables, which, in this case, are elevation changes in the surface due to different processes. 
````R
load('./img/FINAL_ST_review_RL02v15_EEMD/Nuuk_results_0313')

ice_all=subset(x_des,(name=="ICE"))
surf_all=subset(x_des,(name=="SURF_FIRN"))
year = 2003:2013
````

The first step is to create regular meshes to grid our data, as your netcdf product will be easier to create

````R
### start for loop here
dhdt_ice <-matrix(nrow=90000, ncol=10)
dhdt_smb <-matrix(nrow=90000, ncol=10)
std_ice <-matrix(nrow=90000, ncol=10)
std_smb <-matrix(nrow=90000, ncol=10)
for(i in 1:11) {
 i = i-1
number=2003+i

ice=subset(ice_all,(t==i))
surf=subset(surf_all,(t==i))
surf=surf[1:1984,]

G <- interp(ice$x,ice$y,ice$x_mean,xo=seq(-3000,3000,length=300),yo=seq(-3000,3000,length=300))
dhdt_ice[,i]<-G$z

G_std <- interp(ice$x,ice$y,ice$std,xo=seq(-3000,3000,length=300),yo=seq(-3000,3000,length=300))
std_ice[,i]<-G_std$z

GS <- interp(surf$x,surf$y,surf$x_mean,xo=seq(-3000,3000,length=300),yo=seq(-3000,3000,length=300))
dhdt_smb[,i]<-GS$z

GS_std <- interp(surf$x,surf$y,surf$std,xo=seq(-3000,3000,length=300),yo=seq(-3000,3000,length=300))
std_smb[,i]<-GS_std$z
}
````

#image(G,zlim=c(-0.5,0.5))

We use the [raster] (https://cran.r-project.org/web/packages/raster/raster.pdf) package available in R  to create a geospatial dataset and define our projection.

````R
R = raster(G)
projection(R) = CRS("+proj=stere +lat_ts=-71 +datum=WGS84 +units=km")
 ````
Now is time to use the [netcdf](https://cran.r-project.org/web/packages/ncdf4/index.html) package and start packing up our dataset! 
First, we describe our dimensional variables, that is, our X and Y coordinates (or lat and lon, depending on which projection you are working)
````R
Xdim <- ncdim_def("X", "km_polar-stereographic", as.double(G$x))
Ydim <- ncdim_def("Y", "km_polar-stereographic", as.double(G$y))
````
We also define the time slices you have data for
````R
t <- ncdim_def( "Time", "Year", 1:10, unlim=TRUE)
````

Finally, we define our dependent variables
````R
varname="dhdt_ice"
units="m/yr"
dlname <- "dh/dt due to ice dynamic processes"
fillvalue <- -9999
tmp.def <- ncvar_def(varname, units, list(Xdim, Ydim, t), fillvalue, 
                     dlname, prec = "double")
					 
varname="std_ice"
units="m/yr"
dlname <- "std ice dynamics"
fillvalue <- -9999
tmp.def2 <- ncvar_def(varname, units, list(Xdim, Ydim, t), fillvalue, 
                     dlname, prec = "double")
varname="dhdt_smb"
units="m/yr"
dlname <- "dh/dt due to surface processes"
fillvalue <- -9999
tmp.def3 <- ncvar_def(varname, units, list(Xdim, Ydim, t), fillvalue, 
                     dlname, prec = "double")
varname="std_smb"
units="m/yr"
dlname <- "std surface processes"
fillvalue <- -9999
tmp.def4 <- ncvar_def(varname, units, list(Xdim, Ydim, t), fillvalue, 
                     dlname, prec = "double")
````
And let's create the netCDF!
````R
ncfname <- "Final_Results_AISmassbalance_2003-2013.nc"
ncout <- nc_create(ncfname, list(tmp.def,tmp.def2,tmp.def3,tmp.def4))

# put the array
ncvar_put(ncout, tmp.def, dhdt_ice)
ncvar_put(ncout, tmp.def2, std_ice)
ncvar_put(ncout, tmp.def3, dhdt_smb)
ncvar_put(ncout, tmp.def4, std_smb)

# put additional attributes into dimension and data variables
ncatt_put(ncout, "X", "axis", "X")  
ncatt_put(ncout, "Y", "axis", "Y")
 

# add global attributes
title <- "Antarctica annual trends of elevation changes (in m/yr) due surface and ice dynamic processes. Gridded onto a 20 km polar-stereographic"
ncatt_put(ncout, 0, "title", title)

# close the file, writing data to disk
nc_close(ncout)
````