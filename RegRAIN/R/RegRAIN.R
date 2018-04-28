#' RegRAIN Spatial Interpolator for regionalized daily rainfall and air temperature surfaces.
#'
#' This function apply the method RegRAIN to interpolate daily rainfall and air temperature data from climate stations.
#' @param datos Object of class "data.frame" that contains the input rainfall data including the columns:
#'              climate station code,date, longitude, latitude and rainfall data.
#' @param dem  The Digital Elevation Model (DEM) that will be used in the interpolation.
#' @param ini The initial day of interpolation.
#' @param fin The final day of interpolation.
#' @param crossv When "TRUE" a cross validation process will be performed and additional to the RegRAIN raster a plot showing "goodness of fit" indicators will be exported to the working directory. If "FALSE" only the RegRAIN Interpolation raster will be generated.
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
#' LONGITUD:   with the longitudes or X coordinates of the climate station.
#'
#' LATITUD:    with the latitudes or Y coordinates of the climate station.
#'
#' PPT:        with the daily rainfall data in mm (use "." for dec). NA data should be set as "99999" or as "NA" for a proper recognition.
#'
#' Run  View(datos) to see the example data.
#'
#' dem:         is a raster class object. It should contain the projection information and we recommend to use tif format. The extent
#'             of the raster will be the same as the RegRAIN interpolation result extent.
#'
#' ini:         is the initial day of interpolation. Note that the first day is the same as the first date you have in your data object. You should not use dates here, just the number that corresponds
#'             with the initial date you want to interpolate.
#'
#' fin:         is the final day of interpolation. Note that the last day is the same as the final date you have in your data object. You should not use dates here, just the number that corresponds
#'             with the final date you want to interpolate.
#'
#' crossv:      is of class "logical". If set "TRUE" a cross validation process for each RegRAIN raster will be performed. This process will take some aditional time depending on the number of
#'             climate stations.
#'
#' @section About RegRAIN:
#'
#' RegRAIN is a regionalized rain model, that integrates the use of digital elevation models,
#' their transformations (slope and aspect) and homogenized climate information recorded on surface climate stations.
#' The calculation is based on a Multiple Linear Regression (RLM).
#'
#'                         yi = a1*xi1 +a2*xi2 + a3*xi3 + a4*xi4 +a5*xi5 + Ei 	Equation 1
#'
#'    Where,
#'
#'    yi = weather variable in the station "i".
#'
#'    aik = Regression coefficient in the climate station "i" for the weather variable "k".
#'
#'    xik = k variables in the climate station "i" with xi1 = 1.
#'
#'    Ei  = residue.
#'
#'
#'  The 5(k) variables that correspond to geographic and terrain factors are
#'
#'  - Latitude and longitude (decimal degrees or meters).
#'
#'  - Elevation (obtained from DEM in meters).
#'
#'  - Slope and Aspect for each climate station (obtained from DEM in degrees).
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
#' #load the daily rainfall data that will be interpolated
#' datos
#' #Run the RegRAIN function. Here the rainfall data between the
#' #first (ini) and the 5th (fin) date will be interpolated, using
#' #the dem pixel resolution (0.0091 decimal degrees).
#' RegRAIN(datos, dem, 1, 2, crossv = FALSE)
#'
#' @section References:
#'
#' Abteilung Hydrometeorologie. (2013). REGNIE (Regionalisierte Niederschlage): Verfahrensbeschreibung & Nutzeranleitung.
#' Offenbach: Deutscher Wetterdienst - DWD. 9 p.
#'
#' Rauthe, M., Steiner, H., U. Riediger, A., Mazurkiewicz, A., & Gratzki, A. (2013). A Central European precipitacion
#' climatology - Part I: Generation and validation of a high-resolution gridded daily data set (HYRAS). Meteorologische
#' Zeitschrift, 22(3), 235-256.
#'
#' Nunez Lopez, D., Trevino Garza, E., Reyes Gomez, V., Munoz Robles, C., Aguirre Calderon, O., & Jimenez Perez, J. (2014).  
#' Uso de modelos de regresion para interpolar espacialmente la precipitacion media mensual en la cuenca del rio Conchos. 
#' Revista mexicana de ciencias agricolas, 5 (2), 201-213. 
#'
#' @importFrom graphics plot legend lines
#' @importFrom fields Tps
#' @importFrom raster terrain extract extent resample interpolate res<- res writeRaster raster projection crs<- crs
#' @importFrom sp bbox gridded<- gridded proj4string coordinates<- coordinates
#' @importFrom stats complete.cases lm na.omit nls predict residuals
#' @importFrom gstat idw
#' @importFrom grDevices dev.off png
#' @importFrom hydroGOF rmse mae
#' @importFrom utils capture.output
#' @export
#'
#'
RegRAIN <- function(datos,dem,ini,fin,crossv){
  
  if (class(datos$FECHA) !="Date") {#Identificar si la columna fecha en el dataframe con los datos se encuentra en formato "Date"
    
    print(paste("Error: The column FECHA in your dataframe is not of the class Date and in the format %y-%m-%d", sep="")) #Imprime mensaje de error
    
  } else { #Si se encuentra en formato "Date" sigue el proceso
    
    options(digits=20) #Convertir las columnas numericas con 20 digitos
    datos$LONGITUD <- as.numeric(as.character(datos$LONGITUD)) #Convertir la columna de "LONGITUD" al tipo "numeric"
    datos$LATITUD <- as.numeric(as.character(datos$LATITUD)) #Convertir la columna de "LATITUD" al tipo "numeric"
    fechas = unique(datos$FECHA) #Indicar que cada fecha en la columna "FECHA" es un identificador unico
    fechas = sort(fechas) #Ordenar las fechas
    
    #Obtener dataframe con estaciones totales ingresadas
    DATA1 <- datos[!duplicated(datos$CODIGO), ] #Obtener el numero de estaciones a procesar
    coordinates(DATA1) = ~LONGITUD + LATITUD #Convetir a objeto de la clase "SpatialPointsDataFrame"
    
    #DERIVAR ASPECTO DEL DEM
    ASPECTO <- terrain(dem, opt=c('aspect'), unit='degrees', neighbors=8)
    #DERIVAR PENDIENTE DEL DEM
    PENDIENTE <- terrain(dem, opt=c('slope'), unit='degrees', neighbors=8)
    #EXTRAER INFORMACION DE ASPECTO DEL RASTER DE ASPECTO GENERADO
    aspecto <- data.frame(coordinates(DATA1), DATA1$CODIGO, extract(ASPECTO, DATA1)) #Extraer los datos en un data frame
    names(aspecto)[1:4] <- c("X", "Y","CODIGO", "ASP") #Nombrar las columnas
    aspecto <- data.frame(aspecto$CODIGO, aspecto$ASP) #Seleccionar columnas de interes
    names(aspecto)[1:2] <- c("CODIGO", "ASP") #Nombrar las columnas
    aspecto<-aspecto[complete.cases(aspecto),] #Suprimir filas con valores NA del archivo datos
    datos = merge(datos,aspecto,by="CODIGO") #Unir el data frame de extraccion con el archivo de datos totales
    
    #EXTRAER INFORMACION DE ELEVACION DEL RASTER DE ASPECTO GENERADO
    dem1 <- data.frame(coordinates(DATA1), DATA1$CODIGO, extract(dem, DATA1)) #Extraer los datos en un data frame
    names(dem1)[1:4] <- c("X", "Y","CODIGO", "Z") #Nombrar las columnas
    dem1 <- data.frame(dem1$CODIGO, dem1$Z) #Seleccionar columnas de interes
    names(dem1)[1:2] <- c("CODIGO", "DEM") #Nombrar las columnas
    dem1<-dem1[complete.cases(dem1),] #Suprimir filas con valores NA del archivo datos
    datos = merge(datos,dem1,by="CODIGO") #Unir el data frame de extraccion con el archivo de datos totales
    
    #EXTRAER INFORMACION DE PENDIENTES DEL RASTER DE PENDIENTE GENERADO
    pendiente <- data.frame(coordinates(DATA1), DATA1$CODIGO, extract(PENDIENTE, DATA1)) #Extraer los datos en un data frame
    names(pendiente)[1:4] <- c("X", "Y","CODIGO", "PEN") #Nombrar las columnas
    pendiente <- data.frame(pendiente$CODIGO, pendiente$PEN) #Seleccionar columnas de interes
    names(pendiente)[1:2] <- c("CODIGO", "PEN") #Nombrar las columnas
    pendiente<-pendiente[complete.cases(pendiente),] #Suprimir filas con valores NA del archivo datos
    datos = merge(datos,pendiente,by="CODIGO") #Unir el data frame de extraccion con el archivo de datos totales
    
    DATA1<- datos #Crear un data frame de respaldo al de datos totales con las columnas extraidas
    
    p <- raster(extent(bbox(dem))) #Crea un objeto de tipo raster con extent igual al DEM de entrada
    res(p) <- res(dem) #Asigna una resolucion de pixel acorde con el DEM de entrada
    p <- resample(p, dem) #Iguala la estructura y caracteristicas a las del DEM de entrada
    res <- data.frame(res(p)) #Extrae en un data frame la resolucion del DEM de entrada
    res1<-as.numeric(as.character(res[1:1,])) #Extrae como un vector la resolucion en X del DEM
    res2<-as.numeric(as.character(res[2:2,])) #Extrae como un vector la resolucion en Y del DEM
    xy <- data.frame(datos$LATITUD,datos$LONGITUD) #Crea un dataframe con las coordenadas de los datos de entrada
    xy <- unique(xy) #Asigna valores unicos a las coordenadas
    suppressWarnings(spline1 <- Tps(xy,xy$datos.LONGITUD)) #Suprime advertencias y mensajes del proceso de interpolacion
    spline.lat <- interpolate(p,spline1) #Asigna los valores interpolados al raster "p" creado inicialmente en un nuevo raster
    
    suppressWarnings(spline2 <- Tps(xy,xy$datos.LATITUD)) #Suprime advertencias y mensajes del proceso de interpolacion
    spline.long <- interpolate(p,spline2) #Asigna los valores interpolados al raster "p" creado inicialmente en un nuevo raster
    
    #Crear Grid para la interpolacion IDW
    project <- proj4string(dem) #define un objeto de proyeccion geografica a partir del DEM de entrada
    xmin <- xmin(dem) #Extrae xmin del DEM de entrada
    xmax <- xmax(dem) #Extrae xmax del DEM de entrada
    ymin <- ymin(dem) #Exrae ymin del DEM de entrada
    ymax <- ymax(dem) #Extrae ymax del DEM de entrada
    x.range <- as.numeric(c(xmin,xmax)) #Genera un rango de X a partir de xmin y xmax
    y.range <- as.numeric(c(ymin,ymax)) #Genera un rango de Y a partir de ymin y ymax
    grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = res1), y = seq(from = y.range[1], to = y.range[2], by = res2)) #Expande la grilla a partir de los valores minimos y maximos de X - Y
    coordinates(grd) <- ~x + y #Definir coordenadas espaciales para crear un objeto georreferenciado
    gridded(grd) <- TRUE #Asignar tipo grilla
    crs(grd) <- project #Proyectar la grilla
    
    for (i in ini:fin)  { #Bucle de proceso para la interpolacion RegRAIN
      PPT = subset(datos,datos$FECHA==fechas[i]) #Extraccion de la base de datos total para una fecha especifica
      print(paste("Starting RegRAIN Interpolation",fechas[i], sep=" ")) #Imprime mensaje en la consola
      PPT[PPT == 99999] <- NA #Reemplaza todos los valores 99999 por NAs
      PPT = as.data.frame(na.omit(PPT)) #Elimina NAs del analisis
      PPT$FECHA <- as.Date(PPT$FECHA) #Convierte a Date la columna FECHA
      RLM = lm(PPT ~ ASP + DEM + PEN + LONGITUD + LATITUD, data=PPT) #Aplica la Regresion Lineal Multiple
      
      cons <- summary(RLM)$coefficients[1,1] #Extrae la constante de la regresion como vector
      asp  <- summary(RLM)$coefficients[2,1] #Extrae el coeficiente de aspecto como vector
      demc  <- summary(RLM)$coefficients[3,1] #Extrae el coeficiente de elevacion como vector
      pen  <- summary(RLM)$coefficients[4,1] #Extrae el coeficiente de pendiente como vector
      long <- summary(RLM)$coefficients[5,1] #Extrae el coeficiente de longitud como vector
      lat <- summary(RLM)$coefficients[6,1] #Extrae el coeficiente de latitud como vector
      
      #Calcula el residuo de la regresion lineal multiple
      residuo = PPT$PPT - (asp*PPT$ASP) - (demc*PPT$DEM) - (pen*PPT$PEN) - (long*PPT$LONGITUD) - (lat*PPT$LATITUD)
      residuo2 <- as.data.frame(residuo) #convierte el residuo en data frame
      residuo3 <-as.data.frame(na.omit(residuo2)) #Omite NAs
      d <- as.data.frame(row.names(residuo3)) #Prepara el data frame para interpolacion
      residuo4 <- data.frame(residuo3, d) #Prepara el data.frame para interpolacion
      names(residuo4)[1:2] <- c("residuo","ID") #Asigna nombres a columnas
      f <- data.frame(PPT$LONGITUD, PPT$LATITUD, residuo2) #Prepara el data frame para interpolacion
      names(f)[1:3] <- c("X","Y", "residuo") #Asigna nombres a columnas
      coordinates(f) <- ~X + Y #Definir coordenadas espaciales para crear un objeto georreferenciado
      crs(f) <- project #Proyecta el data frame de residuos
      z= capture.output(idw.res <- idw(formula = residuo ~ 1, locations = f, newdata = grd)) # Aplica el interpolador IDW al residuo eliminando mensaje IDW
      raster.idw <- raster(idw.res) #Convierte la superficie IDW a raster
      Residuo <- resample(raster.idw, dem) #Asigna caracteristicas del DEM al raster de IDW generado
      
      RegRAIN_PPT <- (ASPECTO*asp) + (dem*demc) + (PENDIENTE*pen) + (spline.lat*lat) + (spline.long*long) + Residuo #Realiza la interpolacion RegRAIN
      RegRAIN_PPT[RegRAIN_PPT<0] <- 0 #Asigna valores menores a 0 = 0
      
      if (crossv == FALSE) { #Condicional para Cross Validation = FALSE
        
        print(paste("Exporting RegRAIN Raster to working directory", sep="")) #Imprime mensaje en la consola
        writeRaster(RegRAIN_PPT, paste("RegRAIN_PPT_",fechas[i],".tif",sep=""), drivername = "GTiff", overwrite=TRUE) #Exporta raster al directorio de trabajo
        
      } else {
        
        print(paste("Exporting RegRAIN Raster to working directory",sep="")) #Imprime mensaje en la consola
        writeRaster(RegRAIN_PPT, paste("RegRAIN_PPT_",fechas[i],".tif",sep=""), drivername = "GTiff", overwrite=TRUE) #Exporta raster al directorio de trabajo
        
        #Crear data frame para datos de validacion cruzada
        PPT_cv <- PPT
        PPT_cv$PPT_cross <- NA
        PPT_cv[] <- NA
        PPT_cv$FECHA <- as.Date(PPT_cv$FECHA)
        #Crear data frame para datos de validacion cruzada
        codigo = unique(PPT$CODIGO) #Indicar que cada campo en la columna CODIGO es un identificador unico
        codigo = sort(codigo) #Ordenar los codigos
        
        print(paste("Starting Cross Validation Process ",fechas[i],sep="")) #Imprime mensaje en la consola
        
        for (i in 1:nrow(PPT))  {
          CODIGO=NULL #Asigna un identificador a la variable CODIGO
          PPT_cross <- subset(PPT, !(CODIGO %in% codigo[i])) #Extraccion de la base de datos de la fecha interpolada para un codigo especifico
          PPT_cross1 <- subset(PPT, (CODIGO %in% codigo[i])) #Extraccion de la base de datos de la informacion asociada al codigo especifico
          coordinates(PPT_cross1) = ~LONGITUD + LATITUD #Convetir a objeto de la clase "SpatialPointsDataFrame"
          RLM = lm(PPT ~ ASP + DEM + PEN + LONGITUD + LATITUD, data=PPT_cross) #Aplica la Regresion Lineal Multiple
          
          cons <- summary(RLM)$coefficients[1,1] #Extrae la constante de la regresion como vector
          asp  <- summary(RLM)$coefficients[2,1] #Extrae el coeficiente de aspecto como vector
          demc  <- summary(RLM)$coefficients[3,1] #Extrae el coeficiente de elevacion como vector
          pen  <- summary(RLM)$coefficients[4,1] #Extrae el coeficiente de pendiente como vector
          long <- summary(RLM)$coefficients[5,1] #Extrae el coeficiente de longitud como vector
          lat <- summary(RLM)$coefficients[6,1] #Extrae el coeficiente de latitud como vector
          
          #Calcula el residuo de la regresion lineal multiple
          residuo = PPT_cross$PPT - (asp*PPT_cross$ASP) - (demc*PPT_cross$DEM) - (pen*PPT_cross$PEN) - (long*PPT_cross$LONGITUD) - (lat*PPT_cross$LATITUD)
          residuo2 <- as.data.frame(residuo) #Convierte el residuo en data.frame
          residuo3 <-as.data.frame(na.omit(residuo2)) #Omite NAs
          d <- as.data.frame(row.names(residuo3)) #Prepara el data frame para interpolacion
          residuo4 <- data.frame(residuo3, d) #Prepara el data frame para interpolacion
          names(residuo4)[1:2] <- c("residuo","ID") #Asigna nombres a columnas
          f <- data.frame(PPT_cross$LONGITUD, PPT_cross$LATITUD, residuo2) #Prepara el data frame para interpolacion
          names(f)[1:3] <- c("X","Y", "residuo") #Asigna nombres a columnas
          coordinates(f) <- ~X + Y #Definir coordenadas espaciales para crear un objeto georreferenciado
          crs(f) <- project #Proyecta el data frame de residuos
          j= capture.output(idw.res <- idw(formula = residuo ~ 1, locations = f, newdata = grd)) # Aplica el interpolador IDW al residuo eliminando mensaje IDW
          raster.idw <- raster(idw.res) #Convierte la superficie IDW a raster
          Residuo <- resample(raster.idw, dem) #Asigna caracteristicas del DEM al raster de IDW generado
          
          RegRAIN_PPT_cross <- (ASPECTO*asp) + (dem*demc) + (PENDIENTE*pen) + (spline.lat*lat) + (spline.long*long) + Residuo #Realiza la interpolacion RegRAIN
          RegRAIN_PPT_cross[RegRAIN_PPT_cross<0] <- 0 #Asigna valores menores a 0 = 0
          
          PPT_test <- data.frame(coordinates(PPT_cross1), PPT_cross1$CODIGO, extract(RegRAIN_PPT_cross, PPT_cross1)) #Extraer los datos en un data frame
          names(PPT_test)[1:4] <- c("x", "Y","CODIGO", "PPT_cross") #Nombrar las columnas
          PPT_test <- data.frame(PPT_test$CODIGO, PPT_test$PPT_cross) #Seleccionar columnas de interes
          names(PPT_test)[1:2] <- c("CODIGO", "PPT_cross") #Nombrar las columnas
          PPT_test<-PPT_test[complete.cases(PPT_test),]##Suprimir filas con valores NA del archivo datos
          PPT_test1 = merge(PPT,PPT_test,by="CODIGO") #Unir el data frame de extraccion con el archivo de datos totales
          PPT_cv[i,] <- PPT_test1 #Selecciona la fila i del datanfrmae PPT_cv
          fecha <- PPT_cv$FECHA[1] #Extrae la primera fecha del data frame
          avance <- (i/nrow(PPT))*100 #Calcula el avance de la validacion cruzada
          print(paste("Running Cross Validation ",round(avance, digits = 1),"%", sep=" ")) #Imprime mensaje en la consola
          
        }
        fecha <- PPT_cross$FECHA[1] #Extrae la primera fecha del data.frame
        if (sum(PPT_cv$PPT) > 0) { #Condicional que evalua que la suma de la columna PPT sea mayor a "0"
          png(filename= paste("RegRAIN_Cross_Validation_Plot_",fecha,".png",sep=""), #Abre el comando de generacion del png
              units="cm", #Define unidades
              width=10, #Define el ancho del grafico
              height=10, #Define el alto del grafico
              pointsize=6, #Define unidades de tamano de puntos en el grafico
              res=250) #Define resolucion del grafico
          diagram <- plot(PPT_cv$PPT_cross, PPT_cv$PPT, xlab = "Predicted Precipitation (mm)",
                          ylab = "Observed Precipitation(mm)", main= paste("RegRAIN Cross Validation Plot", fecha, sep=" ")) #Crea Plot de Validacion Cruzada
          mod <- nls(PPT_cv$PPT_cross ~ exp(a + b * PPT_cv$PPT), data=PPT_cv, start = list(a = 0, b = 0) ) #Aplica modelo exponencial al grafico
          lines(PPT_cv$PPT, predict(mod, list(x = PPT_cv$PPT),col="red")) # Adiciona curva ajustada al modelo exponencial
          RSS.p<-sum(residuals(mod)^2) #Suma residual de cuadrados
          TSS<-sum((PPT_cv$PPT-mean(PPT_cv$PPT_cross))^2) #Suma total de cuadrados
          r.squared<-1-(RSS.p/TSS) #Calcula coeficiente de determinacion R2
          #Genera leyenda al grafico
          legend("topright", bty="n", legend=c(paste("R2 is", format(r.squared, digits=4)),paste("RMSE is", format(rmse(PPT_cv$PPT_cross, PPT_cv$PPT), digits=4)),paste("MAE is", format(mae(PPT_cv$PPT_cross, PPT_cv$PPT, na.rm=TRUE), digits=4))))
          print(paste("Exporting Cross Validation Plot ",fecha,sep=" ")) #Imprime mensaje en la consola
          dev.off() #Cierra el comando de generacion del png
          print(paste("End of Cross Validation Process ",fecha,sep=" ")) #Imprime mensaje en la consola
        } else {
          png(filename= paste("RegRAIN_Cross_Validation_Plot_",fecha,".png",sep=""), #Abre el comando de generacion del png
              units="cm", #Define unidades
              width=10, #Define el ancho del grafico
              height=10, #Define el alto del grafico
              pointsize=6, #Define unidades de tamano de puntos en el grafico
              res=250) #Define resolucion del grafico
          diagram <- plot(PPT_cv$PPT_cross, PPT_cv$PPT, xlab = "Predicted Precipitation (mm)",
                          ylab = "Observed Precipitation(mm)", main= paste("RegRAIN Cross Validation Plot", fecha, sep=" ")) #Crea Plot de Validacion Cruzada
          #Genera leyenda al grafico
          legend("topright", bty="n", legend=c(paste("RMSE is", format(rmse(PPT_cv$PPT_cross, PPT_cv$PPT), digits=4)),paste("MAE is", format(mae(PPT_cv$PPT_cross, PPT_cv$PPT, na.rm=TRUE), digits=4))))
          dev.off() #Cierra la generacion del png
          print(paste("Exporting Cross Validation Plot ",fecha,sep=" ")) #Imprime mensaje en la consola
          print(paste("End of Cross Validation Process ",fecha,sep=" ")) #Imprime mensaje en la consola
        }
        
      }
      
    }
    
  }
  
}
