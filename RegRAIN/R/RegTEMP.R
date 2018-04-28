#' RegTEMP Spatial Interpolator for regionalized daily air temperature surfaces.
#'
#' This function apply the method RegTEMP to interpolate daily rainfall and air temperature data from climate stations.
#' @param datos Object of class "data.frame" that contains the input temperature data including the columns:
#'              weather station code,date, longitude, latitude anda temperature data.
#' @param dem  The Digital Elevation Model (DEM) that will be used in the interpolation.
#' @param ini The initial day of interpolation.
#' @param fin The final day of interpolation.
#' @param crossv When "TRUE" a cross validation process will be performed and additional to the RegTEMP raster a plot showing "goodness of fit" indicators will be exported to the working directory. If "FALSE" only the RegTEMP Interpolation raster will be generated.
#' @param var Object of class "character" that defines the temperature to interpolate (max, med, min)
#' @details
#' For defining the arguments, the following should be considered:
#'
#'
#' datos       is the object of class "data.frame" containing the daily rainfall data and must
#'             include the following columns:
#'
#' CODIGO:     With the station code.
#'
#' FECHA:      with the date of the data in the format YEAR-MONTH-DAY.
#'             If this column is not of class date and in the format YEAR-MONTH-DAY, the interpolation could be wrong.
#'
#' LONGITUD:   with the longitudes or X coordinates of the weather station.
#'
#' LATITUD:    with the latitudes or Y coordinates of the weather station.
#'
#' TEMP:      with the daily temperature data in degrees Celsius (use "." for dec). NA data should be set as "99999" or as "NA" for a proper recognition.
#'
#' Run  View(datos_tmin) to see the example data.
#'
#' dem:         is a raster class object. It should contain the projection information and we recommend to use *tif format. The extent
#'             of the raster will be the same as the RegTEMP interpolation result extent.
#'
#' ini:         is the initial day of interpolation. Note that the first day is the same as the first date you have in your data object. You should not use dates here, just the number that corresponds
#'             with the initial date you want to interpolate.
#'
#' fin:         is the final day of interpolation. Note that the last day is the same as the final date you have in your data object. You should not use dates here, just the number that corresponds
#'             with the final date you want to interpolate.
#'
#' crossv:      is of class "logical". If set "TRUE" a cross validation process for each RegTEMP raster will be performed. This process will take some aditional time depending on the number of
#'             weather stations.
#'
#'  var         is an object of class "character" that indicates if the temperature is the daily min, max or mean.
#'
#' @section About RegTEMP:
#'
#' RegTEMP is a regionalized air temperature model, that integrates the use of digital elevation models
#  and homogenized climate information recorded on surface climate stations.
#' The calculation of RegTEMP is based on a Multiple Linear Regression (RLM).
#'
#'                         yi = a1*xi1 +a2*xi2 + a3*xi3 + Ei 	Equation 1
#'
#'    Where,
#'
#'    yi = weather variable in the station "i".
#'
#'    aik = Regression coefficient in the weather station "i" for the weather variable "k".
#'
#'    xik = k variables in the weather station "i" with xi1 = 1.
#'
#'    Ei  = residue.
#'
#'
#'  The 3(k) variables that correspond to geographic and terrain factors are
#'
#'  - Latitude and longitude (decimal degrees or meters).
#'
#'  - Elevation (obtained from DEM in meters).
#'
#' @author
#'
#' Diego Alzate.
#' Agroclimatology Group
#' CORPOICA.
#' E-mail: dfalzate@corpoica.org.co
#'
#' @examples
#' #load the digital elevation model to use in the interpolation
#' dem <- raster(system.file("extdata", "demwgs84.tif", package="RegRAIN"))
#' #load the daily minimum tmperature data that will be interpolated
#' datos_tmin
#' #Define wheter it is min, max or mean temperature
#' var<- "Tmin_"
#' #Run the RegTEMP function. Here the minimun temperature data between the
#' #first (ini) and the 5th (fin) date will be interpolated, using
#' #the dem pixel resolution (0.0091 decimal degrees).
#' RegTEMP(datos_tmin, dem, 1, 2, TRUE,"Tmin_")
#'
#' @section References:
#'
#' Abteilung Hydrometeorologie.(2013). REGNIE (Regionalisierte Niederschlaege): Verfahrensbeschreibung & Nutzeranleitung.
#' Offenbach: Deutscher Wetterdienst - DWD. 9 p.
#'
#' Rauthe, M., Steiner, H., U. Riediger, A., Mazurkiewicz, A., & Gratzki, A.(2013). A Central European precipitacion
#' climatology - Part I: Generation and validation of a high-resolution gridded daily data set (HYRAS). Meteorologische
#' Zeitschrift, 22(3), 235-256.
#' 
#' Nunez Lopez, D., Trevino Garza, E., Reyes Gomez, V., Munoz Robles, C., Aguirre Calderon, O., & Jimenez Perez, J. (2014).  
#' Uso de modelos de regresion para interpolar espacialmente la precipitacion media mensual en la cuenca del rio Conchos. 
#' Revista mexicana de ciencias agricolas, 5 (2), 201-213. 
#'
#' @importFrom graphics abline
#'
#' @export
RegTEMP <- function(datos, dem, ini, fin, crossv, var){

  if (class(datos$FECHA) !="Date") {

    print(paste("Error: The column FECHA in your dataframe is not of the class Date and in the format %y-%m-%d", sep=""))

  } else {


  options(digits=20)
  datos$LONGITUD <- as.numeric(as.character(datos$LONGITUD))
  datos$LATITUD <- as.numeric(as.character(datos$LATITUD))
  datos$FECHA <- as.Date(datos$FECHA)
  fechas = unique(datos$FECHA)
  fechas = sort(fechas)

  DATA1 <- datos[!duplicated(datos$CODIGO), ]
  coordinates(DATA1) = ~LONGITUD + LATITUD

  dem1 <- data.frame(coordinates(DATA1), DATA1$CODIGO, extract(dem, DATA1))
  names(dem1)[1:4] <- c("X", "Y","CODIGO", "Z")
  dem1 <- data.frame(dem1$CODIGO, dem1$Z)
  names(dem1)[1:2] <- c("CODIGO", "DEM")
  dem1<-dem1[complete.cases(dem1),]
  datos = merge(datos,dem1,by="CODIGO")

  DATA1<- datos

  p <- raster(extent(bbox(dem)))
  res(p) <- res(dem)
  p <- resample(p, dem)
  res <- data.frame(res(p))
  res1<-as.numeric(as.character(res[1:1,]))
  res2<-as.numeric(as.character(res[2:2,]))
  xy <- data.frame(datos$LATITUD,datos$LONGITUD)
  xy <- unique(xy)
  suppressWarnings(spline1 <- Tps(xy,xy$datos.LONGITUD))
  spline.lat <- interpolate(p,spline1)

  suppressWarnings(spline2 <- Tps(xy,xy$datos.LATITUD))
  spline.long <- interpolate(p,spline2)

  project <- proj4string(dem)
  xmin <- xmin(dem)
  xmax <- xmax(dem)
  ymin <- ymin(dem)
  ymax <- ymax(dem)
  x.range <- as.numeric(c(xmin,xmax))
  y.range <- as.numeric(c(ymin,ymax))

  grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = res1), y = seq(from = y.range[1], to = y.range[2], by = res2))

  #Definir coordenadas espaciales para crear un objeto georreferenciado
  coordinates(grd) <- ~x + y
  gridded(grd) <- TRUE
  crs(grd) <- project

  for (i in ini:fin)  {
    TEMP = subset(datos, datos$FECHA==fechas[i])
    print(paste("Starting RegTEMP Interpolation",fechas[i], sep=" "))
    TEMP[TEMP == 99999] <- NA
    TEMP = as.data.frame(na.omit(TEMP))
    TEMP$FECHA <- as.Date(TEMP$FECHA)
    RLM = lm(TEMP ~ DEM + LONGITUD + LATITUD, data=TEMP)

    cons <- summary(RLM)$coefficients[1,1]
    demc  <- summary(RLM)$coefficients[2,1]
    long <- summary(RLM)$coefficients[3,1]
    lat <- summary(RLM)$coefficients[4,1]


    residuo = TEMP$TEMP - (demc*TEMP$DEM) - (long*TEMP$LONGITUD) - (lat*TEMP$LATITUD)
    residuo2 <- as.data.frame(residuo)
    residuo3 <-as.data.frame(na.omit(residuo2))
    d <- as.data.frame(row.names(residuo3))
    residuo4 <- data.frame(residuo3, d)
    names(residuo4)[1:2] <- c("residuo","ID")
    f <- data.frame(TEMP$LONGITUD, TEMP$LATITUD, residuo2)
    names(f)[1:3] <- c("X","Y", "residuo")
    coordinates(f) <- ~X + Y
    crs(f) <- project
    z= capture.output(idw.res <- idw(formula = residuo ~ 1, locations = f, newdata = grd))
    raster.idw <- raster(idw.res)
    Residuo <- resample(raster.idw, dem)

    RegTEMP_TEMP <- (dem*demc) + (spline.lat*lat) + (spline.long*long) + Residuo

    if (crossv == FALSE) {

      print(paste("Exporting RegTEMP Raster to working directory", sep=""))
      writeRaster(RegTEMP_TEMP, paste("RegTEMP_",var, fechas[i],".tif",sep=""), drivername = "GTiff", overwrite=TRUE)

    } else {

      print(paste("Exporting RegTEMP Raster to working directory",sep=""))
      writeRaster(RegTEMP_TEMP, paste("RegTEMP_",var, fechas[i],".tif",sep=""), drivername = "GTiff", overwrite=TRUE)

      #Crear data.frame para datos de validacion cruzada
      TEMP_cv <- TEMP
      TEMP_cv$TEMP_cross <- NA
      TEMP_cv[] <- NA
      TEMP_cv$FECHA <- as.Date(TEMP_cv$FECHA)
      #Crear data.frame para datos de validacion cruzada
      codigo = unique(TEMP$CODIGO)
      codigo = sort(codigo)

      print(paste("Starting Cross Validation Process ",fechas[i],sep=""))

      for (i in 1:nrow(TEMP))  {
        CODIGO=NULL
        TEMP_cross <- subset(TEMP, !(CODIGO %in% codigo[i]))
        TEMP_cross1 <- subset(TEMP, (CODIGO %in% codigo[i]))
        coordinates(TEMP_cross1) = ~LONGITUD + LATITUD
        RLM = lm(TEMP ~ DEM + LONGITUD + LATITUD, data=TEMP_cross)

        cons <- summary(RLM)$coefficients[1,1]
        demc  <- summary(RLM)$coefficients[2,1]
        long <- summary(RLM)$coefficients[3,1]
        lat <- summary(RLM)$coefficients[4,1]

        residuo = TEMP_cross$TEMP - (demc*TEMP_cross$DEM) - (long*TEMP_cross$LONGITUD) - (lat*TEMP_cross$LATITUD)
        residuo2 <- as.data.frame(residuo)
        residuo3 <-as.data.frame(na.omit(residuo2))
        d <- as.data.frame(row.names(residuo3))
        residuo4 <- data.frame(residuo3, d)
        names(residuo4)[1:2] <- c("residuo","ID")
        f <- data.frame(TEMP_cross$LONGITUD, TEMP_cross$LATITUD, residuo2)
        names(f)[1:3] <- c("X","Y", "residuo")
        #Definir coordenadas espaciales para crear un objeto georreferenciado
        coordinates(f) <- ~X + Y
        crs(f) <- project
        j= capture.output(idw.res <- idw(formula = residuo ~ 1, locations = f, newdata = grd))
        raster.idw <- raster(idw.res)
        Residuo <- resample(raster.idw, dem)


        RegTEMP_TEMP_cross <- (dem*demc) + (spline.lat*lat) + (spline.long*long) + Residuo

        TEMP_test <- data.frame(coordinates(TEMP_cross1), TEMP_cross1$CODIGO, extract(RegTEMP_TEMP_cross, TEMP_cross1))
        names(TEMP_test)[1:4] <- c("x", "Y","CODIGO", "TEMP_cross")
        TEMP_test <- data.frame(TEMP_test$CODIGO, TEMP_test$TEMP_cross)
        names(TEMP_test)[1:2] <- c("CODIGO", "TEMP_cross")
        TEMP_test<-TEMP_test[complete.cases(TEMP_test),]##Suprimir filas con valores NA del archivo datos
        TEMP_test1 = merge(TEMP,TEMP_test,by="CODIGO")
        TEMP_cv[i,] <- TEMP_test1
        fecha <- TEMP_cv$FECHA[1]
        avance <- (i/nrow(TEMP))*100
        print(paste("Running Cross Validation ",round(avance, digits = 1), "%", sep=" "))
      }

      fecha <- TEMP_cross$FECHA[1]
      if (sum(TEMP_cv$TEMP) > 0) {
        png(filename= paste("RegTEMP_Cross_Validation_Plot_",fecha,".png",sep=""),
            units="cm",
            width=10,
            height=10,
            pointsize=6,
            res=250)
        diagram <- plot(TEMP_cv$TEMP_cross, TEMP_cv$TEMP, xlab = paste("Predicted Temperature","(",var,")","(Celsius)", sep=" "),
                        ylab = paste("Observed Temperature","(",var,")","(Celsius)",sep=" "), main= paste("RegTEMP Cross Validation Plot", fecha, sep=" "))
        abline(fit <- lm(TEMP_cv$TEMP_cross ~ TEMP_cv$TEMP, data=TEMP_cv), col='red')
        legend("bottomright", bty="n", legend=c(paste("R2 is", format(summary(fit)$adj.r.squared, digits=4)),paste("RMSE is", format(rmse(TEMP_cv$TEMP_cross, TEMP_cv$TEMP), digits=4)),paste("MAE is", format(mae(TEMP_cv$TEMP_cross, TEMP_cv$TEMP, na.rm=TRUE), digits=4))))
        dev.off()
        print(paste("Exporting Cross Validation Plot ",fecha,sep=" "))
        print(paste("End of Cross Validation Process ",fecha,sep=" "))
      } else {
        png(filename= paste("RegTEMP_Cross_Validation_Plot_",fecha,".png",sep=""),
            units="cm",
            width=10,
            height=10,
            pointsize=6,
            res=250)
        diagram <- plot(TEMP_cv$TEMP_cross, TEMP_cv$TEMP, xlab = paste("Predicted Temperature","(",var,")","(Celsius)", sep=" "),
                        ylab = paste("Observed Temperature","(",var,")","(Celsius)",sep=" "), main= paste("RegTEMP Cross Validation Plot", fecha, sep=" "))
        legend("bottomright", bty="n", legend=c(paste("RMSE is", format(rmse(TEMP_cv$TEMP_cross, TEMP_cv$TEMP), digits=4)),paste("MAE is", format(mae(TEMP_cv$TEMP_cross, TEMP_cv$TEMP, na.rm=TRUE), digits=4))))
        dev.off()
        print(paste("Exporting Cross Validation Plot ",fecha,sep=" "))
        print(paste("End of Cross Validation Process ",fecha,sep=" "))
       }

     }

   }

 }

}


