test<- full_simulation_archive[[12]]

coordinates(test) <- ~x + y

plot(test)

x_range<- as.numeric(range(test$x))
y_range<- as.numeric(range(test$y))

grd <- expand.grid(x = seq(from = x_range[1],
                           to = x_range[2], 
                           by = 0.1),
                   y = seq(from = y_range[1],to = y_range[2], 
                           by = 0.1))  # expand points to grid
# Convert grd object to a matrix and then turn into a spatial
# points object
coordinates(grd) <- ~x + y
# turn into a spatial pixels object
gridded(grd) <- TRUE

#### view grid with points overlayed
plot(grd, cex = 1.5, col = "grey")
plot(test,
     pch = 15,
     col = "red",
     cex = 1,
     add = TRUE)

# interpolate the data
idw_pow1 <- idw(formula = age_hours ~ 1,
                locations = test,
                newdata = grd,
                idp = 1)

test_ras<- rast(idw_pow1)

plot(test_ras,
     col = matlab.like(30))

library(terra)

# 1. Convert points to a SpatVector
pts <- vect(test, geom=c("x", "y"), crs="EPSG:4326")

# 2. Create a buffer around each point (set width to your desired distance)
pts_buffer <- buffer(pts, width = 0.2) 

# 3. Aggregate individual overlapping buffers into a single combined mask
mask_poly <- aggregate(pts_buffer)

# 4. Mask the interpolated raster
masked_raster <- mask(x=test_ras,mask =  mask_poly)

# 5. Plot the result
plot(masked_raster$var1.pred,col=matlab.like(30))
points(test$x,test$y,col="darkgreen",pch=16)



