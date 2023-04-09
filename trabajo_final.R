#####  Librerías a utilizar ####
library('geoR')
library('spdep')
library('gstat')
library("leafsync")
library('mapview')
library('leaflet')
library('RColorBrewer')
library('ggplot2')
library('ggmap')
library('tibble')
library('caret')
library('sf')
library('sp')
library('PerformanceAnalytics')
library("rgdal")
library('lattice') 
library('grDevices')
library('GGally')
library('dplyr')

#####  Carga de información ####
data(meuse)
str(meuse)

copia_seguridad <- meuse

#### Hacemos una breve exploración de la base ####
# Vemos la cantidad de registros y columnas
dim(meuse[,c('x','y','zinc')]) # nos quedan 155 observaciones con la longitud, latitud y valor de Zinc reportado

# Vemos algunos datos de la base
head(meuse[,c('x','y','zinc')])
View(meuse[,c('x','y','zinc')])

# Vemos la estadística descriptiva
summary(meuse[,c('x','y','zinc')]) # Tenemos valores máximos que superan ampliamente la mediana de la población. 
# Sospechamos entonces que habrán eventuales datos atípicos

# Analizamos duplicados y valores nulos
sum(duplicated(meuse[,c('x','y','zinc')])) # No tenemos registros duplicados
sum(is.na(meuse[,c('x','y','zinc')])) # No tenemos valores NA. 

# OBS: Son importantes estos chequeos porque en caso de haber duplicados o nulos tendríamos que
# evaluar si imputar o eliminar registros.

# Definimos una función para hacer el formato del gráfico, que tendrá incluída}
# una regresión lineal clásica y una regresión generada por observaciones locales.

my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p
}
# Graficamos correlaciones y nube de puntos
g = ggpairs(meuse[,c('x','y','zinc')],columns = 1:3, lower = list(continuous = my_fn))
g

# Vemos que parece ser que el valor de Zinc disminuye para mayores niveles de longitud reportada (es
# decir, para mayores valores de x), mientras que parece crecer a mayores niveles de latitud
# reportada (es decir, a mayores niveles de y). Sin embargo, dichas tendencias parecen ser muy
# leves, y con un gran nivel de varianza a lo largo de sus rangos. Las correlaciones en efecto
# tienen valores próximos a 0 y sin significancia estadística alguna.

# Analizamos la normalidad de los datos
x <- meuse$x
y <- meuse$y
zinc<-meuse$zinc

df = as.data.frame(zinc)

ggplot(df, aes(x=zinc)) + geom_histogram() + labs(title = "Distribución del Zinc",y = 'Frequency') + theme(plot.title = element_text(hjust=0.5)) 
#Podemos observar que los datos no siguen una distribución Normal, más bien parecerían tener 
#una distribución Exponencial.

ggplot(df, aes(x=zinc)) + geom_boxplot() + labs(title = "Boxplot del Zinc") + theme(plot.title = element_text(hjust=0.5)) 
#Volvemos a ver la asimetría y valores que parecen ser atípicos.

ggplot(df, aes(sample=zinc)) + stat_qq() + stat_qq_line() + labs(title = 'QQ Plot de la distribución del Zinc') + theme(plot.title = element_text(hjust=0.5))
#Vemos que se aleja de los cuantiles de una Normal

# Test de Normalidad
shapiro.test(meuse[,c('zinc')]) # Tenemos información suficiente como para rechazar la hipótesis
# nula de que los datos siguen una distribución normal.

# Vemos como parece quedar con una distribución más cercanaa a la normal si lo llevamos a logaritmo
meuse$lnZn <- log(meuse$zinc, base = exp(1))
hist(meuse$lnZn, breaks = 16)
shapiro.test(meuse$lnZn)
# Igual rechazamos normalidad

# Hacemos un plot básico:
# Creamos un Data Frame espacial
coordinates(meuse) <- c("x", "y")

plot(meuse["zinc"], asp = 1, pch = 1, cex = 4*meuse$zinc/max(meuse$zinc), col=heat.colors(3))
# Agregamos los márgenes del río
data(meuse.riv)
lines(meuse.riv)
# Vemos las mayores concentraciones cerca del río

# De manera más elegante, usamos "MapView". Para esto, volvemos a subir la base ya que necesitamos
# reformatear el CRS
data(meuse)
coordinates(meuse) <- ~x+y
proj4string(meuse) <- CRS("+init=epsg:28992")

mapview(meuse, zcol = c("zinc"), legend = TRUE)

# Podemos ver múltiples cosas

m1 <- mapview(meuse, zcol = "zinc", legend = TRUE)
m2 <- mapview(meuse, zcol = "landuse", map.types = "Esri.WorldImagery", legend = TRUE)
m3 <- mapview(meuse, zcol = "dist.m", legend = TRUE)
m4 <- mapview(meuse, zcol = "elev", legend = TRUE)
sync(m1, m2, m3, m4) # 4 paneles sincronizados

x<-meuse$x
y<-meuse$y
zinc<-meuse$zinc
data<-cbind(x,y,zinc)
datag<-as.geodata(data)
plot(datag)

#A través de este análisis exploratorio podemos entender que los datos no tienen una distribución normal, sino más
#bien algo exponencial.
#Podemos observar que todos los puntos se encuentran dispersos a lo largo del río, teniendo ciertos lugares donde 
#hay más concentración de los mismos, como es el caso de los puntos azules. Por otro lado, en el caso de los puntos rojos,
#notamos que se encuentran dispersos, casi que en su totalidad, en el lado izquierdo del río, donde anterirormente encontramos
#una mayor contaminación de zinc. 
#En el caso de los puntos verdes, notamos que se encuentrar dispersos por el medio del río, donde referenciandolo con 
#la localización previamente realizada podemos notar que serían los puntos de menor contaminación por zinc.
#Por último, identificando los puntos amarillos, encontramos que estos son los que hay en menor medida,denotando estos un grado
#menor de contaminación por zinc. 
#Vemos un río ancho y no tan largo, porque este río es conocido como el "más" ancho de Europa
#Haciendo referencia a la relación entre la variable en estudio con cada coordenada, podemos notar una tendencia en la relación
#de la contaminación por zinc con respecto a la coordenada x e y, pero de una magnitud muy chica y sin significancia
#estadística. 

#######################################################
### Índice de Moran Global, Local e Índice de Geary ###
######### Estadísticos de Moran y Geary ###############
#######################################################

# Volvemos a cargar los datos porque tenemos que volver a formatearle las coordenadas
#Vamos a generar una grilla para ver definir el vecindario y ver los vecinos de cada punto

data(meuse)
coordinates(meuse) <- c("x", "y")
coordenadas <- coordinates(meuse)

# Recolectamos las coordenadas x - y de todos los puntos. 
grilla <- dnearneigh(meuse,0,400)
card(grilla)
plot(grilla, coordenadas, pch=18, cex=0.5)

# Asignamos pesos espaciales
pesos <- nb2listw(grilla, style = "W")

# Ploteamos el vecindario
plot(grilla, coordenadas, col = "red", pch = 19, cex = 1)

#Calculamos el los índices de Moran local y global, y el índice de Geary

### Índice de Moran ###

#Con el índice de Moran vamos a querer evaluar la autocorrelación espacial para poder comprender la variación de un fénomeno,
#en este caso del zinc en el río Meuse, en un marco geográfico de análisis.
#Para decir que hay evidencia de exitencia de autocorrelación positiva deberíamos notar que el zinc se agrupa en zonas uniformes
#conformando de esta manera una especie de cluster. En el caso de la autocorrelación negativa, deberíamos notar que el zinc
#se encuentra disperso, es decir que es la presencia de zinc será disímil en sus logares aledaños/vecinos. 
#Por último, una posible conclusión podría ser que no encontramos autocorrelación espacial, esta puede ser arribada cuando
#encontramos un comportamiento aleatorio en el estudio del fénomeno. 

### Índice Moran Global###

#Generamos un gráfico que evalúa cuan similar es cada dato con respecto a los datos de sus vecinos 
#A través de este índice podemos encontrar datos atípicos, es decir posibles outliers espaciales. 

# Generamos un gráfico que evalue cuan similar es el dato respecto a sus vecinos
M <- moran.plot(meuse$zinc,pesos,zero.policy=F,col=3, quiet=TRUE,labels=T,xlab = "zinc", ylab="lag(zinc)")
View(M)
# Tenemos un montón de potenciales outliers, como por ejemplo la observación 118 y la 69.

if (require(ggplot2, quietly=TRUE)) {
  xname <- attr(M, "xname")
  ggplot(M, aes(x=x, y=wx)) + geom_point(shape=1) + 
    geom_smooth(formula=y ~ x, method="lm") + 
    geom_hline(yintercept=mean(M$wx), lty=2) + 
    geom_vline(xintercept=mean(M$x), lty=2) + theme_minimal() + 
    geom_point(data=M[M$is_inf,], aes(x=x, y=wx), shape=9) +
    geom_text(data=M[M$is_inf,], aes(x=x, y=wx, label=labels, vjust=1.5)) +
    xlab('Zinc') + ylab(paste0("Spatially lagged ", 'Zinc'))
}

### Índice de Moran Local ###
#A través del índice de Moran Local vamos a poder detectar la presencia de posibles outliers espaciales, ayudandonos a elimiar
#de la muestra los datos que difieren significativamente de su vecindario, en otras palabras cuatificará la similitud
#del punto con respecto a la vecindad.
#Este índice será calculado para cada observación.

# Calculamos el índice de Moran Local y mostramos los resultados
ML <- localmoran(meuse$zinc, pesos, alternative ="less")
IML <-printCoefmat(data.frame(ML,row.names=elevation$Casos),check.names=FALSE)

# Potencial outlier, puesto que por su valor de autocorrelación espacial, se trata de una
# observación que tiene un valor que no se correlaciona mucho con su entorno.
IML[69,]

# Acá vemos el listado ordenado en función de su IML
head(IML[order(IML$Pr.z...E.Ii..),])

# Generamos una matriz que muestra cualitativamente 
# los datos en relación a sus vecinos:
# FALSE: dato similar a su entorno.
Elim = M$is_inf
Elim
#Como podemos ver las observaciones que son iguales a TRUE son las que demarcamos como posibles outliers

### Test de Moran ####
#Con el test de Moran queremos determinar la presencia o ausencia significativa
#de autocorrelación espacial dentro de un conjunto de datos.
#Para este caso la hipótesis nula planteara la no correlación espacial de los datos, por lo que si el p valor obtenido
#del test es menor que el nivel de significancia elegido, entonces la hipótesis nula será rechazada, concluyendo
#la presencia de autocorrelación espacial en los datos.

# Hacemos el test de Moran Global (IMG) de asociación espacial
moran.test(meuse$zinc, nb2listw(grilla, style = "W"))
# Hay pruebas suficientes para rechazar la hipótesis nula de que "No hay correlación
# espacial entre las variables que generan los datos"

moran.test(meuse$zinc, nb2listw(grilla, style = "S"))
moran.test(meuse$zinc, nb2listw(grilla, style = "B"))
moran.test(meuse$zinc, nb2listw(grilla, style = "C"))
moran.test(meuse$zinc, nb2listw(grilla, style = "U"))
moran.test(meuse$zinc, nb2listw(grilla, style = "minmax"))

# Rechazamos la hipótesis nula tanto para los pesos que usan una codificación binaria ("B"),
# como para los que usan pesos estandarizados ("W"), para los que usan pesos globalmente
# estandarizados ("C"), para los que usan los pesos globalmente estandarizados divididos
# por el tamaño de la vecindad, para los que tienen una codificación que permite estabilizar
# la varianza ("S") y para los que estandarizan los pesos con el rango de máximo y mínimo
# ("minmax").

### Índice de Geary ###

#Con el índice de Geary vamos a medir la asociación espacial a través del cambio, es decir vamos a analizar la asociación
#espacial a través de la varianza del cambio.
#A través de este índice evaluaremos la distribución del zinc en el río Meuse es de manera aleatoria o no.
#Diremos que hay una alta autocorrelación espacial si la diferencia de los valores del zinc en diferentes lugares son
#similares y la distancia entre los lugares es pequeña.
#Este índice es más sensible a outliers que el índice de Moran.
#El índice de Geary se expresa como un valor entre 0 y 2, donde 1 indica una distribución espacial aleatoria 
#de la variable. Un valor menor que 1 indica una autocorrelación positiva, lo que significa que los valores similares 
#tienden a estar cerca unos de otros, mientras que un valor mayor que 1 indica una autocorrelación negativa, 
#lo que significa que los valores diferentes tienden a estar cerca unos de otros.

# Hacemos el teste de Geary Local (Varianza del cambio)
geary.test(meuse$zinc, nb2listw(grilla, style = "W"),randomisation = FALSE)
geary.test(meuse$zinc, nb2listw(grilla, style = "S"),randomisation = FALSE)
geary.test(meuse$zinc, nb2listw(grilla, style = "B"),randomisation = FALSE)
geary.test(meuse$zinc, nb2listw(grilla, style = "C"),randomisation = FALSE)
geary.test(meuse$zinc, nb2listw(grilla, style = "U"),randomisation = FALSE)
geary.test(meuse$zinc, nb2listw(grilla, style = "minmax"),randomisation = FALSE) 

# Con todos las pruebas de pesos, el valor del estadístico es menor a 1, 
# lo cual nos indica que en sitios conectados los valores del zinc son similares.

# ¿Qué pasa si hacemos el test de Moran sacando Outliers?
moran.test(meuse$zinc, nb2listw(grilla, style = "W")) # Con todo

# Prueba sin potenciales outliers: Obervación 57, 125, 69, 122 y 98
sin_outlier = copia_seguridad[-c(57,125,69,122,98),]

coordinates(sin_outlier) <- c("x", "y")
coordenadas_sin_outlier <- coordinates(sin_outlier)

# Recolectamos las coordenadas x - y de todos los puntos. 
grilla_sin_outlier <- dnearneigh(sin_outlier,0,400)
card(grilla_sin_outlier)
plot(grilla_sin_outlier, coordenadas_sin_outlier, pch=18, cex=0.5)
moran.test(sin_outlier$zinc, nb2listw(grilla_sin_outlier, style = "W")) 
geary.test(sin_outlier$zinc, nb2listw(grilla_sin_outlier, style = "W"),randomisation = FALSE)
# Vemos que ahora ambos índices se mueven a favor de una hipótesis de mayor relación espacial,
# aunque no es un cambio espectacular.
#Al arrojar valores menores a 1, podemos afirmar que hay autocorrelación positiva, es decir que los valores similares
#tienden a estar cerca unos de otros. 


#### Variogramas ####
# Esto es lo que respecta al análisis estructural.

#El variograma es una herramienta para el análisis geoestadístico utilizada para describir la variabilidad espacial
#de un fenómeno, en este caso el zinc. 


# Genero el Variograma Nube
#A través del variograma nube vamos a poder graficar la distancia de los puntos en el espacio vs la variable de cambio (semivarianza al cuadrado)
nube_clasica <- variog(meuse, coords = coordinates(meuse), data = meuse$zinc, option = "cloud")
nube_CH <- variog(meuse, coords = coordinates(meuse), data = meuse$zinc, option = "cloud", estimator.type = "modulus")

# Genero el variograma experimental/ empírico
#Como con el variograma nube no alcanza vamos a calcular el variograma empírico, el cual va a ser construido a partir del variograma nube, dividiendo el eje de distancias en 
#intervalos, tomando un representante sobre el cual promedio con todos los puntos que caen dentro del intervalo, consiguiendo así la representación del mismo.
bin_clasico <- variog(meuse, coords = coordinates(meuse),uvec = seq(0, 1000, by = 50), data = meuse$zinc)
bin_CH <- variog(meuse, coords = coordinates(meuse), data = meuse$zinc, uvec = seq(0, 1000, by = 50), estimator.type= "modulus")

# Lo ploteamos
par(mfrow=c(2,2))
plot(nube_clasica, main = "classical estimator")
plot(nube_CH, main = "modulus estimator")
plot(bin_clasico, main = "classical estimator")
plot(bin_CH, main = "modulus estimator")
par(mfrow = c(1,1))

#A través de estos gráficos podríamos observar un ligero comportamiento anisotropíco en la varble, dado que para algunos puntos el comportamiento del zinc
#presenta algunas variaciones en ciertas direcciones, pero las mismas son muy pequeñas. 


### Boxplots de bins

bin1 <- variog(meuse, coords = coordinates(meuse), data = meuse$zinc, uvec = seq(0, 1000, by = 50), bin.cloud = T)
bin2 <- variog(meuse, coords = coordinates(meuse), data = meuse$zinc, uvec = seq(0, 1000, by = 50), estimator.type = "modulus", bin.cloud = T)

par(mfrow = c(1,2))
plot(bin1, bin.cloud = T, main = "classical estimator")
plot(bin2, bin.cloud = T, main = "modulus estimator")
par(mfrow = c(1,2))


## Datos direccionales

#Calculamos el variograma variando las direcciones, de esta manera queremos entender la variabilidad espacial del zinc en diferentes direcciones. 
#En este caso queremos entender si estamos frente a una vairable anisotrópica o isotrópica, es decir queremos terminar de confirmar si el comportamiento
#estadístico varía según la dirección en la que se mida la distancia. 
#Para poder calcularlo vamos a tomar un ángula y distancia determinados

vario.2 <- variog(meuse, coords = coordinates(meuse), data = meuse$zinc, uvec = seq(0, 1000, by = 50), dir=0)
vario.3 <- variog(meuse, coords = coordinates(meuse), data = meuse$zinc, uvec = seq(0, 1000, by = 50), dir=pi/2)
vario.4 <- variog(meuse, coords = coordinates(meuse), data = meuse$zinc, uvec = seq(0, 1000, by = 50), dir=pi/4)

par(mfrow = c(1,1))
plot(bin_clasico, type="l")
lines(vario.2, lty = 2, col = 2)
lines(vario.3, lty = 3, col = 3)
lines(vario.4, lty = 4, col = 4)
legend("bottomright", c("omnidireccional", "0°", "90°", "45°"), col=c(1,2,3,4), lty=c(1,2,3,4))


# Otro plot
varias_direcciones = variog4(coords = coordinates(meuse), data = meuse$zinc, uvec = seq(0, 1000, by = 50))
plot(varias_direcciones)

#Podemos observar que para diferentes distancias y ángulos la variable en cuestión presenta ciertas variaciones, por lo que estaríamos en presencia de anisotrópia. 


#Variogramas de residuos
res1.v = variog(meuse, coords = coordinates(meuse), data = meuse$zinc, uvec = seq(0, 1000, by = 50))
set.seed(123)
plot(res1.v)

# Intervalos de simulación por permutación aleatoria de los residuos
s1 = variog.mc.env(meuse, coords = coordinates(meuse), data = meuse$zinc, obj = res1.v)
plot(res1.v, env = s1)

# Usando una simulación simple de Monte Carlo, vemos que en ambos variograsmas empíricos (osea, tanto el de los residuos como
# con el clásico sin tendencia) nos muestras que hay puntos que caen fuera del "envelope", lo cual nos dice que hay margen para estudiar la depen-
# dencia espacial, sobre todo a distancias cortas.

s2 = variog.mc.env(meuse, coords = coordinates(meuse), data = meuse$zinc, obj = bin_clasico)
plot(bin_clasico, env = s2)

## AJUSTE de un  Variograma a los datos
# A ojo

res1.v.ef = eyefit(res1.v)
sigmasq = 143151.79
phi = 333.33
nugget = 17893.97
modelo = "exp"

plot(res1.v)

# Reemplazamos curva, ajustando el mejor Variograma Teórico que se ajuste
lines.variomodel(res1.v,
                 cov.model = modelo, 
                 cov.par= c(sigmasq,phi),
                 nug = nugget)

# Ahora vamos a calcular el variograma sin tendencia con la función "variogram"
v <- variogram(zinc~1, meuse) # Sin tendencia, pues no hayamos una
plot(v)
v

# Con estos valores podemos dar unos valores iniciales para que el modelo estime
# el variograma teórico
#Calculamos el variograma teórico porque necesitamos de una línea continúa para poder estimar las observaciones que no fueron representadas. 

vt_exp = fit.variogram(v, vgm(190000, "Exp", 1400, 30000))
vt_exp # Nugget de 9486, el rango es 9486.599 + 163285.464, aunque después fija
# el rango en 381.7098 (el cuál tendremos que multiplicar por 3)
rango_practico = 381.7098 * 3
rango_practico
plot(v , vt_exp)

# Probamos un variograma esférico
vt_sph = fit.variogram(v, vgm(190000, "Sph", 1400, 30000))
vt_sph # Acá al rango no habría que considerarlo
plot(v , vt_sph)

# ¿Cúal ajusta mejor? El de menor sumas del cuadrado del error
attr(vt_exp, 'SSErr')
attr(vt_sph, 'SSErr')
# Con este criterio el que mejor ajusta es el primero.

# Chequeamos como nos queda suponiendo una tendencia lineal
vten <- variogram(zinc~x+y, meuse)
plot(vten)

# Ajustamos a este variograma un modelo teórico exponencial, con tendencia
vtent_exp = fit.variogram(vten, vgm(190000, "Exp", 1400, 30000))
vtent_exp
plot(vten , vtent_exp)

# Ahora ajustamos un modelo te?rico esf?rico, con tendencia
vtent_sph = fit.variogram(vten, vgm(190000, "Sph", 1400, 30000))
vtent_sph
plot(vten , vtent_sph)

# Comparación de los modelos con tendencia propuestos 
attr(vtent_exp, 'SSErr') # Sigue ajustando mejor el exponencial, cómo en el caso previo. 
attr(vtent_sph, 'SSErr')

# Ahora bien, ¿vamos con o sin tendencia? No podemos mezclar los error cuadráticos
# medio entre modelos con y sin tendencia. Tenemos que ver los plots, pero como vimos 
# antes, no hay patrones claros.
g
# Además, al haber hecho variogramas con y sin tendencia y obtener resultados tan parecidos,
# se refuerza nuestra hipótesis de que no tienen tendencia.

# ¿Son isotrópicos?
# Cutoff es la maxima distancia que tiene sentido para considerar un distancia entre los puntos. Consideramos 4200.

v <- variogram(zinc~1, meuse, cutoff = 4200, width = 50, map=T)
plot(v)
# Parece un proceso Anisotrópico, con dirección por continuidad, alcanzando otros
# valores si corremos la dirección perpendicular al sentido de la mancha
v1 <- variogram(zinc~1, meuse, cutoff = 4200, width = 200, map=T)
plot(v1)
# Con un mayor ancho de ventana parece más evidente

show.vgms()


vExp <- fit.variogram(vten, vgm(model = "Exp"),fit.method = 2)
#vGau <- fit.variogram(vten, vgm(model = "Gau"),fit.method = 2)
vSph <- fit.variogram(vten, vgm(model = "Sph"),fit.method = 2)
vMat <- fit.variogram(vten, vgm(model = "Mat", nugget = 1,kappa = 0.5),fit.method = 2)
vBes <- fit.variogram(vten,vgm("Bes"),fit.method = 2)
vSte <- fit.variogram(vten,vgm("Ste"),fit.method = 2)

d <- c("Exponencial"= plot ( vten , vExp),
       #"Gaussiano"=plot ( vten , vGau),
       "Esferico "=plot ( vten , vSph),
       "Matern"=plot ( vten , vMat),
       "Stein's"= plot ( vten , vSte),
       "Bessel"= plot ( vten , vBes))
d

vExpLine=variogramLine(vExp,500)
#vGauLine=variogramLine(vGau,500)
vSphLine=variogramLine(vSph,500)
vMatLine=variogramLine(vMat,500)
vSteLine=variogramLine(vSte,500)
vBesLine=variogramLine(vBes,500)

ggplot(mapping = aes(dist,gamma))+
  geom_point(data = vten)+
  geom_line(data = vExpLine,aes(color="Exponencial"))+
  geom_line(data = vSphLine,aes(color="Esferico"))+
  geom_line(data = vMatLine,aes(color="Matern"))+
  geom_line(data = vSteLine,aes(color="Stein's"))+
  geom_line(data = vBesLine,aes(color="Bessel"))+
  scale_color_discrete("Modelo")+
  theme_classic()


#### Kriging ####
#El objetivo del Kriging es generar una predicción basado en el variograma teórico para el punto donde se desconoce el valor de la variable,
#es decir queremos generar una predicción suave donde a través de estas podamos predecir el valor de lo que se desconoce, es decir el valor de la variable donde no fue medida.
#El Kriging nos dará una predicción lineal, la cual estima el valor desconocido a través de los referentes más cercanos a la variable de interés. 
#Decimos que es lineal porque es la suma de los valores en las posiciones vecinas, po lo que es una correlación lineal, donde se le agrega peso en base a la influencia que tienen
#sobre la variable. 


# Kriging Ordinario: Se aplica para Procesos estacionarios con media desconocida #
# AJUSTE CON VARIOGRAMA EMPÍRICO SIN TENDENCIA #

# Cómo la librería trabaja con un objeto geoespacial distinto a los casos previos, volvemos a tomar la base
datos <- select(copia_seguridad,x,y,zinc)
coordinates(datos) = ~x+y
datosg<-as.geodata(datos) # Convertimos en datos geospaciales sin usar GeoR, sino como GeoData

# Calculamos y graficamos un semivariograma
vv <- variogram(zinc~1, meuse, cutoff = 4200, width = 50, map=T)
plot(vv)
vv <- variogram(zinc~1, meuse, cutoff = 4200, width = 500)
plot(vv)
vv <- variogram(zinc~1, meuse,cloud=T) # Sin tendencia, pues no hayamos una
plot(vv)


# Para ajustar el modelo, volvemos a hacer un variograma empírico
v <- variogram(zinc~1, meuse) # Sin tendencia, pues no hayamos una
plot(v)

# Ajustamos varios modelos de variograma
# y seleccionamos dos

# Ajuste 1 (Exp)
v1_exp = fit.variogram(v, vgm(190000, "Exp", 1400, 30000))
v1_exp
plot(v , v1_exp) #Podemos ver un efecto sill o meceta en el valor 700, esto quiere decir que a partir de este valor obtenemosd el valor de variación máxima.
#Por otro lado podemos observar un efecto nugget o pepita. 

# Ajuste 2 (Sph)

v2_sph = fit.variogram(v, vgm(190000, "Sph", 1400, 30000))
v2_sph
plot(v, v2_sph)

# Ajuste 3 (Gau). Tenemos que estimar los parámetros de inicio primero
res1.v.gaus = eyefit(res1.v)
res1.v.gaus
v3_gau = fit.variogram(v, vgm(143151, "Gau", 333.3, 4))
v3_gau
plot(v , v3_gau)


attr(v1_exp, 'SSErr')
attr(v2_sph, 'SSErr')
attr(v3_gau, 'SSErr')

# Nos quedamos con el modelo exponencial, que es el que menor error nos arroja.

# Sin embargo, también probamos un cuarto ajuste automático con varios modelos.

# Ajuste 4
v4_sel = fit.variogram(v, vgm(143151, c("Exp", "Sph", "Mat"), 333.3, 4))
v4_sel
plot(v , v4_sel)
attr(v4_sel, 'SSErr')

#Vamos a quedarnos con los dos mejores modelos
# Sph y Exp (seleccionado por R)

-------------------------------------------------
#Generamos una grilla de predicción llamada g2
#Elegimos posiciones en el espacio donde no tenemos datos
points(datosg)
g2 <- expand.grid(x=seq(179500, 180500, by=100), y=seq(330000,332000, by=100))
gridded(g2) = ~x+y
plot(g2)

# Realizamos la predicción sobre la grilla g2 con
# los dos modelos competidores: 

# Exponencial
k1 <- krige(zinc~1, datos, g2, model = v1_exp, nmax = 155)
# Esférico
k2 <- krige(zinc~1, datos, g2, model = v2_sph, nmax = 155)

# Ploteamos el kriging
# Modelo exponencial
spplot(k1["var1.pred"], main = "Kriging ordinario: Valores Predichos (Exp)", col.regions=terrain.colors)
spplot(k1["var1.var"],  main = "Kriging ordinario: Varianza de las Predicciones (Exp)", col.regions=terrain.colors)

# Modelo esférico
spplot(k2["var1.pred"], main = "Kriging ordinario: Valores Predichos (Sph)", col.regions=terrain.colors)
spplot(k2["var1.var"],  main = "Kriging ordinario: Varianza de las Predicciones (Sph)", col.regions=terrain.colors)

# Vemos que la varianza en las predicciones del modelo exponencial parece ser más baja respecto al segundo modelo propuesto.

# Tabla_1 (Modelo Exponencial)
Predicciones1 = k1$var1.pred
Varianza1 = k1$var1.var

# Armamos la tabla
x1=k1$x
y1=k1$y
Tabla_1=cbind(x1,y1, Predicciones1, Varianza1)


# Tabla_2 (Modelo Esférico)
Predicciones2 = k2$var1.pred
Varianza2 = k2$var1.var

# Armamos la tabla
x2=k2$x
y2=k2$y
Tabla_2 = cbind(x2,y2, Predicciones2, Varianza2)

# Hacemos una validación cruzada de los modelos
# Recordamos los dos modelos de variograma ajustados
# y sus parámetros

modelo1 <- vgm(190000, "Exp", 1400, 30000) # Modelo Exponencial
modelo2 <- vgm(134743, "Sph", 24798, 4.73) # Modelo Esférico

# Hacemos validación cruzada del modelo 1
valcruz1 <- krige.cv(zinc~1, datos, modelo1, nfold=155)
valcruz1
names(valcruz1)

# Hacemos validación cruzada del modelo 2
valcruz2 <- krige.cv(zinc~1, datos, modelo2, nfold=155)
valcruz2
names(valcruz2)

# Error medio de predicción.
# Se espera que sea lo mas proximo a cero posible.
mean(valcruz1$residual)
mean(valcruz2$residual)
# En error medio de predicción, el modelo exponencial sigue siendo el gran ganador.

# Error cuadratico medio de predicción.
# Valors bajos son mejores.
mean(valcruz1$residual^2)
mean(valcruz2$residual^2)
# En ECM de predicción, el modelo esférico le gana al exponencial.

# Error cuadrático medio normalizado.
# Valor deseado: proximo a 1.
mean(valcruz1$zscore^2)
mean(valcruz2$zscore^2)
# En el error cuadrático medio normalizado el modelo exponencial gana. 

# Correlación lineal entre valores observados y predichos
cor(valcruz1$observed, valcruz1$observed - valcruz1$residual)
cor(valcruz2$observed, valcruz2$observed - valcruz2$residual)

# Ploteamos dicha correlación
par(mfrow=c(1,2))
plot(valcruz1$observed,valcruz1$observed - valcruz1$residual,xlab="Observados (Exp)", ylab="Predichos (Exp)")
plot(valcruz2$observed,valcruz2$observed - valcruz2$residual,xlab="Observados (Sph)", ylab="Predichos (Sph)")
par(mfrow=c(1,1))

# La correlación lineal entre los valores observados y los predichos toma valores bastante similares en 
# ambos modelos.

# Salidas de regresión
# Modelo 1: Exponencial
r1<-valcruz1$observed - valcruz1$residual
regresion1 <- lm(valcruz1$observed ~ r1, data = valcruz1)
summary(regresion1)

# Modelo 2: Esférico
r2<-valcruz2$observed - valcruz2$residual
regresion2 <- lm(valcruz2$observed ~ r2, data = valcruz2)
summary(regresion2)

# Tanto el R2 ajustado como los coeficientes obtenidos son muy similares entre ambos modelos, aunque nosotros
# estamos más a favor del caso 1 (exponencial) por la evidencia previa. Los coeficientes estimados en ambas 
# regresiones son similares, siendo el resultante del modelo esférico con un R2 ajustado mínimamente superior.

# Kriging Simple: Se aplica para Procesos estacionarios con media conocida #
# AJUSTE CON VARIOGRAMA EMPÍRICO SIN TENDENCIA #

# Tomamos los datos que venimos usando. Para el caso de la grilla, en este caso tomaremos
# la grilla que viene con la base original. 

# Cargamos la grilla
data(meuse.grid)
gridded(meuse.grid) = ~x+y

# Realizamos la predicción sobre la nueva grilla con los dos modelos competidores: 
# Previamente, tenemos que imputar la media del proceso en el terreno, ya que de otra
# manera estaríamos haciendo un kriging oridnario.

media_zinc = mean(datos$zinc)

# Modelo Exponencial
k1_ord <- krige(zinc~1, datos, meuse.grid, model = v1_exp, nmax = 155, beta = media_zinc)
# Modelo Esférico
k2_ord <- krige(zinc~1, datos, meuse.grid, model = v2_sph, nmax = 155, beta = media_zinc)

# Ahora realizamos las predicciones sobre la grilla nueva para ambos modelos: 

# Modelo Exponencial
spplot(k1_ord["var1.pred"], main = "Kriging simple: Valores Predichos (Exp)", col.regions=terrain.colors)
spplot(k1_ord["var1.var"],  main = "Kriging simple: Varianza de las Predicciones (Exp)", col.regions=terrain.colors)

# Modelo Esférico
spplot(k2_ord["var1.pred"], main = "Kriging simple: Valores Predichos (Sph)", col.regions=terrain.colors)
spplot(k2_ord["var1.var"],  main = "Kriging simple: Varianza de las Predicciones (Sph)", col.regions=terrain.colors)

# Vemos que ambos modelos ajustan bastante bien, ya que tenemos los mayores valores predichos
# muy cercano a lo que es la superficie del río. A su vez, la varianza de las predicciones se
# ven bastante similares en ambos casos.

# Tabla_1_ord (Modelo Exponencial)
Predicciones1_ord = k1_ord$var1.pred
Varianza1_ord = k1_ord$var1.var

# Armamos la tabla
x1_ord=k1_ord$x
y1_ord=k1_ord$y
Tabla_1_ord=cbind(x1_ord,y1_ord, Predicciones1_ord, Varianza1_ord)

# Tabla_2_ord (Modelo Esférico)
Predicciones2_ord = k2_ord$var1.pred
Varianza2_ord = k2_ord$var1.var

# Armamos la tabla
x2_ord=k2_ord$x
y2_ord=k2_ord$y
Tabla_2_ord = cbind(x2_ord,y2_ord, Predicciones2_ord, Varianza2_ord)

# Hacemos una validación cruzada de los modelos
# Hacemos validación cruzada del modelo 1
valcruz1_ord <- krige.cv(zinc~1, datos, modelo1, nfold=155)

# Hacemos validación cruzada del modelo 2
valcruz2_ord <- krige.cv(zinc~1, datos, modelo2, nfold=155)

# Error medio de predicción.
# Se espera que sea lo mas proximo a cero posible.
mean(valcruz1_ord$residual)
mean(valcruz2_ord$residual)
# En error medio de predicción, el modelo exponencial sigue siendo el gran ganador.

# Error cuadratico medio de predicción.
# Valors bajos son mejores.
mean(valcruz1_ord$residual^2)
mean(valcruz2_ord$residual^2)
# En ECM de predicción, el modelo esférico le gana al exponencial.

# Error cuadrático medio normalizado.
# Valor deseado: proximo a 1.
mean(valcruz1_ord$zscore^2)
mean(valcruz2_ord$zscore^2)
# En el error cuadrático medio normalizado el modelo exponencial gana. 

# Correlación lineal entre valores observados y predichos
cor(valcruz1_ord$observed, valcruz1_ord$observed - valcruz1_ord$residual)
cor(valcruz2_ord$observed, valcruz2_ord$observed - valcruz2_ord$residual)

# Ploteamos dicha correlación
par(mfrow=c(1,2))
plot(valcruz1_ord$observed,valcruz1_ord$observed - valcruz1_ord$residual,xlab="Observados (Exp)", ylab="Predichos (Exp)")
plot(valcruz2_ord$observed,valcruz2_ord$observed - valcruz2_ord$residual,xlab="Observados (Sph)", ylab="Predichos (Sph)")
par(mfrow=c(1,1))

# La correlación lineal entre los valores observados y los predichos toma valores muy similares,
# siendo más similares los resultados que en el caso del kriging ordinario.

# Salidas de regresión
# Modelo 1: Exponencial
r1_ord <-valcruz1_ord$observed - valcruz1_ord$residual
regresion1_ord <- lm(valcruz1_ord$observed ~ r1_ord, data = valcruz1_ord)
summary(regresion1_ord)

# Modelo 2: Esférico
r2_ord<-valcruz2_ord$observed - valcruz2_ord$residual
regresion2_ord <- lm(valcruz2_ord$observed ~ r2_ord, data = valcruz2_ord)
summary(regresion2_ord)

# Los resultados son muy similares en comparación con el caso del kriging ordinario. Tanto
# los R2 ajustados como los coeficientes estimados y los p-value son en ambos casos muy parecidos.

# No consideramos necesario hacer Kriging Universal puesto que no encontramos ninguna tendencia
# entre el zinc y alguna de las coordenadas.

# Kriging Universal: Se aplica para Procesos estacionarios con media desconocida, suponiendo que hay una tendencia del zinc con las coordenadas #
# AJUSTE CON VARIOGRAMA EMPÍRICO CON TENDENCIA #
# Vamos a hacerlo de cualquier manera, para evidenciar que los resultados son similares, y fortalecer nuestra hipótesis de que no existe
# una tendencia.

# Tomamos los datos que venimos usando. Volvemos a tomar la grilla que viene con la base original. 

# Cargamos la grilla
data(meuse.grid)
gridded(meuse.grid) = ~x+y

# Realizamos la predicción sobre la nueva grilla con los dos modelos competidores: 

# Modelo Exponencial
k1_uni <- krige(zinc~1, datos, meuse.grid, model = v1_exp, nmax = 155)
# Modelo Esférico
k2_uni <- krige(zinc~1, datos, meuse.grid, model = v2_sph, nmax = 155)

# Ahora realizamos las predicciones sobre la grilla nueva para ambos modelos: 

# Modelo Exponencial
spplot(k1_uni["var1.pred"], main = "Kriging simple: Valores Predichos (Exp)", col.regions=terrain.colors)
spplot(k1_uni["var1.var"],  main = "Kriging simple: Varianza de las Predicciones (Exp)", col.regions=terrain.colors)

# Modelo Esférico
spplot(k2_uni["var1.pred"], main = "Kriging simple: Valores Predichos (Sph)", col.regions=terrain.colors)
spplot(k2_uni["var1.var"],  main = "Kriging simple: Varianza de las Predicciones (Sph)", col.regions=terrain.colors)

# Vemos que ambos modelos ajustan bastante bien, ya que tenemos los mayores valores predichos
# muy cercano a lo que es la superficie del río. Sin embargo, vemos que a nivel general es muy similar al modelo ordinario sin 
# tendencia, lo cual robustece nuestra postura inicial de que no encontramos una tendencia clara entre el zinc y las coordenadas.

# Tabla_1_uni (Modelo Exponencial)
Predicciones1_uni = k1_uni$var1.pred
Varianza1_uni = k1_uni$var1.var

# Armamos la tabla
x1_uni=k1_uni$x
y1_uni=k1_uni$y
Tabla_1_uni=cbind(x1_uni,y1_uni, Predicciones1_uni, Varianza1_uni)

# Tabla_2_uni (Modelo Esférico)
Predicciones2_uni = k2_uni$var1.pred
Varianza2_uni = k2_uni$var1.var

# Armamos la tabla
x2_uni=k2_uni$x
y2_uni=k2_uni$y
Tabla_2_uni = cbind(x2_uni,y2_uni, Predicciones2_uni, Varianza2_uni)

# Hacemos una validación cruzada de los modelos
# Hacemos validación cruzada del modelo 1
valcruz1_uni <- krige.cv(zinc~1, datos, modelo1, nfold=155)

# Hacemos validación cruzada del modelo 2
valcruz2_uni <- krige.cv(zinc~1, datos, modelo2, nfold=155)

# Error medio de predicción.
# Se espera que sea lo mas proximo a cero posible.
mean(valcruz1_uni$residual)
mean(valcruz2_uni$residual)
# En error medio de predicción, el modelo exponencial sigue siendo el gran ganador.

# Error cuadratico medio de predicción.
# Valores bajos son mejores.
mean(valcruz1_uni$residual^2)
mean(valcruz2_uni$residual^2)
# En ECM de predicción, el modelo esférico le gana al exponencial.

# Error cuadrático medio normalizado.
# Valor deseado: proximo a 1.
mean(valcruz1_uni$zscore^2)
mean(valcruz2_uni$zscore^2)
# En el error cuadrático medio normalizado el modelo exponencial gana. 

# Correlación lineal entre valores observados y predichos
cor(valcruz1_uni$observed, valcruz1_uni$observed - valcruz1_uni$residual)
cor(valcruz2_uni$observed, valcruz2_uni$observed - valcruz2_uni$residual)

# Ploteamos dicha correlación
par(mfrow=c(1,2))
plot(valcruz1_uni$observed,valcruz1_uni$observed - valcruz1_uni$residual,xlab="Observados (Exp)", ylab="Predichos (Exp)")
plot(valcruz2_uni$observed,valcruz2_uni$observed - valcruz2_uni$residual,xlab="Observados (Sph)", ylab="Predichos (Sph)")
par(mfrow=c(1,1))

# Obtenemos resultados muy similares a los que obtuvimos previamente con el cómputo del kriging ordinario.

# Salidas de regresión
# Modelo 1: Exponencial
r1_uni <-valcruz1_uni$observed - valcruz1_uni$residual
regresion1_uni <- lm(valcruz1_uni$observed ~ r1_uni, data = valcruz1_uni)
summary(regresion1_uni)

# Modelo 2: Esférico
r2_uni<-valcruz2_uni$observed - valcruz2_uni$residual
regresion2_uni <- lm(valcruz2_uni$observed ~ r2_uni, data = valcruz2_uni)
summary(regresion2_uni)

#### Material Útil ####
# Páginas
# https://rpubs.com/daniballari/geostat_basico
# https://rpubs.com/quarcs-lab/spatial-autocorrelation
# http://www.leg.ufpr.br/lib/exe/fetch.php/patrick:spatialcourse:slides.pdf

# Video. 1:50:00 de la segunda parte de la clase 4 habla un poco de cómo encarar un análisis de estadística espacial.