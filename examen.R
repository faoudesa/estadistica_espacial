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
hist(meuse[,c('zinc')], breaks = 16) # No parecen bastante normales.

boxplot(meuse[,c('zinc')]) # Volvemos a ver la asimetría y valores que parecen ser atípicos.

qqnorm(meuse[,c('zinc')]) # Se aleja de los cuantiles de una normal.
qqline(meuse[,c('zinc')], col=2)

# Test de Normalidad
shapiro.test(meuse[,c('zinc')]) # Tenemos información suficiente como para rechazar la hipótesis
# nula de que los datos siguen una distribución normal.

# Vemos como queda bastante normal si lo llevamos a logaritmo
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

#### Índices Locales y Globales ####
# Volvemos a cargar los datos porque tenemos que volver a formatearle las coordenadas

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

# Generamos un gráfico que evalue cuan similar es el dato respecto a sus vecinos
M <- moran.plot(meuse$zinc,pesos,zero.policy=F,col=3, quiet=T,labels=T,xlab = "zinc", ylab="lag(zinc)")
View(M)
# Tenemos un montón de potenciales outliers, como por ejemplo la observación 118 y la 69.

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

# Hacemos el teste de Geary Local
geary.test(meuse$zinc, nb2listw(grilla, style = "W"),randomisation = FALSE)
geary.test(meuse$zinc, nb2listw(grilla, style = "S"),randomisation = FALSE)
geary.test(meuse$zinc, nb2listw(grilla, style = "B"),randomisation = FALSE)
geary.test(meuse$zinc, nb2listw(grilla, style = "C"),randomisation = FALSE)
geary.test(meuse$zinc, nb2listw(grilla, style = "U"),randomisation = FALSE)
geary.test(meuse$zinc, nb2listw(grilla, style = "minmax"),randomisation = FALSE) 
# Con todos las pruebas de pesos, el valor del estadístico es menor a 1, 
# lo cual nos indica que en sitios conectados los valores del zinc son similares.

# Qué pasa si hacemos el test de Moran sacando Outliers?
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

#### Variogramas ####
# Esto es lo que respecta al análisis estructural.

# Genero el Variograma Nube
nube_clasica <- variog(meuse, coords = coordinates(meuse), data = meuse$zinc, option = "cloud")
nube_CH <- variog(meuse, coords = coordinates(meuse), data = meuse$zinc, option = "cloud", estimator.type = "modulus")

# Genero el variograma experimental
bin_clasico <- variog(meuse, coords = coordinates(meuse),uvec = seq(0, 1000, by = 50), data = meuse$zinc)
bin_CH <- variog(meuse, coords = coordinates(meuse), data = meuse$zinc, uvec = seq(0, 1000, by = 50), estimator.type= "modulus")

# Lo ploteamos
par(mfrow=c(2,2))
plot(nube_clasica, main = "classical estimator")
plot(nube_CH, main = "modulus estimator")
plot(bin_clasico, main = "classical estimator")
plot(bin_CH, main = "modulus estimator")
par(mfrow = c(1,1))

### Boxplots de bins

bin1 <- variog(meuse, coords = coordinates(meuse), data = meuse$zinc, uvec = seq(0, 1000, by = 50), bin.cloud = T)
bin2 <- variog(meuse, coords = coordinates(meuse), data = meuse$zinc, uvec = seq(0, 1000, by = 50), estimator.type = "modulus", bin.cloud = T)

par(mfrow = c(1,2))
plot(bin1, bin.cloud = T, main = "classical estimator")
plot(bin2, bin.cloud = T, main = "modulus estimator")
par(mfrow = c(1,2))

## Datos direccionales

vario.2 <- variog(meuse, coords = coordinates(meuse), data = meuse$zinc, uvec = seq(0, 1000, by = 50), dir=0)
vario.3 <- variog(meuse, coords = coordinates(meuse), data = meuse$zinc, uvec = seq(0, 1000, by = 50), dir=pi/2)
vario.4 <- variog(meuse, coords = coordinates(meuse), data = meuse$zinc, uvec = seq(0, 1000, by = 50), dir=pi/4)

par(mfrow = c(1,1))
plot(bin_clasico, type="l")
lines(vario.2, lty = 2, col = 2)
lines(vario.3, lty = 3, col = 3)
lines(vario.4, lty = 4, col = 4)
legend("bottomright", c("omnidireccional", "0", "90", "45"), col=c(1,2,3,4), lty=c(1,2,3,4))

# Otro plot
varias_direcciones = variog4(coords = coordinates(meuse), data = meuse$zinc, uvec = seq(0, 1000, by = 50))
plot(varias_direcciones)

#Variogramas de residuos
res1.v = variog(meuse, coords = coordinates(meuse), data = meuse$zinc, uvec = seq(0, 1000, by = 50))
set.seed(123)
plot(res1.v)

# Intervalos de simulación por permutación aleatoria de los residuos
s1 = variog.mc.env(meuse, coords = coordinates(meuse), data = meuse$zinc, obj = res1.v)
plot(res1.v, env = s1)

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
# Parece un proceso Anisotrópico, parece un dirección por continuidad, alcanzando otros
# valores si corremos la dirección perpendicular al sentido de la mancha
v1 <- variogram(zinc~1, meuse, cutoff = 4200, width = 200, map=T)
plot(v1)
# Con un mayor ancho de ventana parece más evidente
  

#### Material Útil ####
# Páginas
# https://rpubs.com/daniballari/geostat_basico
# https://rpubs.com/quarcs-lab/spatial-autocorrelation