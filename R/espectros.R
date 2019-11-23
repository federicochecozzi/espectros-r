library(tidyverse)
library(pheatmap)
library(gplots)
library(pracma)
source("R/cargaarchivos.R")
getwd()

.nombre.archivos <- dir(path = "data/", pattern = "[PRT]")
.nombre.archivos
datos.muestras <- cargar.muestras(.nombre.archivos, "data/")

#añado columna con el grupo al que pertenece cada muestra (un grupo corresponde a un compuesto)
datos.muestras <- mutate(datos.muestras, grupo = factor(str_remove(archivo,"[123]")))
glimpse(datos.muestras)#puede observarse las columnas redundantes/con nombres arbitrarios

#limpieza de columnas, la tabla de datos posee muchas columnas redundantes (se puede apreciar en los warnings)
datos.muestras <- datos.muestras %>% 
  select(`Wavelength [nm]`,Intensity,intensity,intensity_1,intensity_2,intensity_3,intensity_4,intensity_5,intensity_6,intensity_7,intensity_8,intensity_9,archivo,grupo) %>% 
  rename(Measure1 = Intensity,Measure2 = intensity,Measure3 = intensity_1,Measure4 = intensity_2,Measure5 = intensity_3,Measure6 = intensity_4,Measure7 = intensity_5,Measure8 = intensity_6,Measure9 = intensity_7,Measure10 = intensity_8,Measure11 = intensity_9, Wavelength = `Wavelength [nm]`)


#preguntas iniciales: ¿necesitamos varios espectros por muestra? ¿qué correlación tienen los espectros?
#las variables relevantes son las longitudes de onda, el número de medición/barrido es un identificador
datos.tidy <-  datos.muestras %>% 
  gather(key = "Nbarrido", value = "Intensidad", Measure1, Measure2, Measure3,Measure4, Measure5, Measure6,Measure7, Measure8, Measure9,Measure10, Measure11) %>%
  spread(key = Wavelength, value = "Intensidad") %>%
  select(grupo,everything()) %>% #reordenando las variables para mayor comodidad
  #luego de eso remuevo los barridos que experimentalmente se consideraron de baja calidad
  filter(((archivo == "PE1" | archivo == "PE2") & (Nbarrido == "Measure3" | Nbarrido == "Measure4" | Nbarrido == "Measure5") ) | ((archivo != "PE1" & archivo != "PE2") & (Nbarrido == "Measure2" | Nbarrido == "Measure3" | Nbarrido == "Measure4")))
  
#etiquetando para facilitar trabajo futuro
rownames(datos.tidy)<-str_c(datos.tidy$archivo,"-B",str_remove(datos.tidy$Nbarrido,"Measure"))

#primero debería comparar espectros entre distintas señales para observar posibles anomalías
#x e y serían longitud de onda e intesidad, así que necesito ordenar los datos un poco en un nuevo tibble
datos.gathered <- datos.tidy %>% gather(key = "Wavelength", value = "Intensity",-grupo,-archivo,-Nbarrido, convert = TRUE)
datos.gathered <- datos.gathered %>% mutate(barrido = str_c(archivo,"-B",str_remove(Nbarrido,"Measure")))
#la idea es que pueda visualizarse los espectros para compararlos, compuestos similares deberían tener el mismo espectro
datos.gathered %>% filter(grupo == "PE") %>% arrange(Wavelength) %>%
  ggplot(aes(x=Wavelength,y=Intensity,color=archivo )) + 
  geom_line(aes(alpha = 0.1))

datos.gathered %>% filter(grupo == "Pu") %>%
  ggplot(aes(x=Wavelength,y=Intensity,color=archivo )) + 
  geom_line(aes(alpha = 0.1))

datos.gathered %>% filter(grupo == "RD") %>%
  ggplot(aes(x=Wavelength,y=Intensity,color=archivo )) + 
  geom_line(aes(alpha = 0.1))

datos.gathered %>% filter(grupo == "TN") %>%
  ggplot(aes(x=Wavelength,y=Intensity,color=archivo )) + 
  geom_line(aes(alpha = 0.1))

#en otros análisis me interesa tratar cada longitud de onda como una variable y que cada barrido sea un vector fila
#calculo la correlación entre diferentes barridos, como cada barrido ocupa una fila tengo que transponer
#nótese que la información en el tibble está ordenada por grupo y archivo ya
scor <- datos.tidy %>% select_if(is.numeric) %>% t() %>% cor()

pheatmap(scor) #heatmap de las correlaciones; demasiadas observaciones para ser cómodo

#es más fácil de visualizar si comparo entre dos grupos
scorPEPu <- datos.tidy %>% rownames_to_column('barrido') %>%
  filter(grupo == "PE" | grupo == "Pu") %>% column_to_rownames('barrido') %>% 
  select_if(is.numeric) %>% t() %>% cor()

pheatmap(scorPEPu)


scorPERD <- datos.tidy %>% rownames_to_column('barrido') %>%
  filter(grupo == "PE" | grupo == "RD") %>% column_to_rownames('barrido') %>% 
  select_if(is.numeric) %>% t() %>% cor()

pheatmap(scorPERD)

scorPETN <- datos.tidy %>% rownames_to_column('barrido') %>%
  filter(grupo == "PE" | grupo == "TN") %>% column_to_rownames('barrido') %>% 
  select_if(is.numeric) %>% t() %>% cor()

pheatmap(scorPETN)

scorPuRD <- datos.tidy %>% rownames_to_column('barrido') %>%
  filter(grupo == "Pu" | grupo == "RD") %>% column_to_rownames('barrido') %>% 
  select_if(is.numeric) %>% t() %>% cor()

pheatmap(scorPuRD)

scorPuTN <- datos.tidy %>% rownames_to_column('barrido') %>%
  filter(grupo == "Pu" | grupo == "TN") %>% column_to_rownames('barrido') %>% 
  select_if(is.numeric) %>% t() %>% cor()

pheatmap(scorPuTN)

scorRDTN <- datos.tidy %>% rownames_to_column('barrido') %>%
  filter(grupo == "RD" | grupo == "TN") %>% column_to_rownames('barrido') %>% 
  select_if(is.numeric) %>% t() %>% cor()

pheatmap(scorRDTN)

#PCA
df_pca <- datos.tidy %>% select(-grupo,-archivo,-Nbarrido) %>% prcomp(scale. = T)
df_out <- as.data.frame(df_pca$x)
labels <- datos.tidy %>% select(grupo,archivo,Nbarrido) %>% rownames_to_column('barrido') 
df_out <- df_out %>% rownames_to_column('barrido') %>% 
  left_join(labels, by = 'barrido') %>% column_to_rownames('barrido')

screeplot(df_pca)#las primeras dos componentes son suficientes

ggplot(df_out,aes(x=PC1,y=PC2,color=grupo )) + geom_point() + xlim(-20,20)


ggplot(df_out,aes(x=PC1,y=PC2,color=Nbarrido )) + geom_point()

#extracción de características útiles para el análisis espectral
#recomendado por compañeros de trabajo
f = findintervals#busca el índice del número más cercano al provisto en un vector numérico
#lo utilizo para obtener el índice de una longitud de onda dada
featurevector <- datos.gathered %>% group_by(barrido) %>% 
  mutate(val = Intensity[f(481.8,Wavelength)], Intensity = Intensity/val) %>%
  summarise(w1 = f(471,Wavelength),w2 = f(474.5,Wavelength),
            integral471a474.5 = trapz(Wavelength[w1:w2],Intensity[w1:w2]),
            w1 = f(488,Wavelength),w2 = f(497,Wavelength),
            integral488a497 = trapz(Wavelength[w1:w2],Intensity[w1:w2]),
            w1 = f(498,Wavelength),w2 = f(506,Wavelength),
            integral498a506 = trapz(Wavelength[w1:w2],Intensity[w1:w2]),
            w1 = f(509,Wavelength),w2 = f(515,Wavelength),
            integral509a515 = trapz(Wavelength[w1:w2],Intensity[w1:w2]),
            pico472.6 = Intensity[f(472.6,Wavelength)],
            pico492.99 = Intensity[f(492.99,Wavelength)],
            pico500.7 = Intensity[f(500.7,Wavelength)],
            pico512.24 = Intensity[f(512.24,Wavelength)]
            ) %>%
  select(-w1,-w2) %>% 
  column_to_rownames('barrido')

featurepca <- prcomp(featurevector,scale. = TRUE)
featureout <- as.data.frame(featurepca$x)
featureout <- featureout %>% rownames_to_column('barrido') %>% 
  left_join(labels, by = 'barrido') %>% column_to_rownames('barrido')

screeplot(featurepca)#las primeras dos componentes son suficientes

ggplot(featureout,aes(x=PC1,y=PC2,color=grupo )) + geom_point() + xlim(-1.5,1.5)
