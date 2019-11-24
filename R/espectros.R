library(pheatmap)
library(gplots)
library(pracma)
library(gridExtra)
library(tidyverse)
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
  mutate(barrido = str_c(archivo,"-B",str_remove(Nbarrido,"Measure"))) %>% 
  group_by(barrido) %>% 
  filter(between(Wavelength,465,520)) %>% #longitudes de onda de interés, determinada quien me suministró los espectros (¿será una buena elección?)
  mutate(
         val = interp1(Wavelength,Intensidad,481.8),#val es temporal
         Intensidad = Intensidad/val #normalización
         ) %>% 
  spread(key = Wavelength, value = "Intensidad") %>%
  select(-val) %>%
  #select(grupo,everything()) %>% #reordenando las variables para mayor comodidad
  #luego de eso remuevo los barridos que experimentalmente se consideraron de baja calidad
  filter((archivo %in% c("PE1", "PE2") & Nbarrido %in% c("Measure3", "Measure4", "Measure5")) | (!(archivo %in% c("PE1", "PE2")) & Nbarrido %in% c("Measure2", "Measure3", "Measure4")))
  
#etiquetando para facilitar trabajo futuro
datos.tidy <- column_to_rownames(datos.tidy, 'barrido')

#primero debería comparar espectros entre distintas señales para observar posibles anomalías
#x e y serían longitud de onda e intesidad, así que necesito ordenar los datos un poco en un nuevo tibble
datos.gathered <- datos.tidy %>% 
  gather(key = "Wavelength", value = "Intensity",-grupo,-archivo,-Nbarrido, convert = TRUE)

datos.gathered <- datos.gathered %>% 
  mutate(barrido = str_c(archivo,"-B",str_remove(Nbarrido,"Measure")))

#la idea es que pueda visualizarse los espectros para compararlos, compuestos similares deberían tener el mismo espectro
p1 <- datos.gathered %>% filter(grupo == "PE") %>% arrange(Wavelength) %>%
  ggplot(aes(x=Wavelength,y=Intensity,color=archivo )) + 
  geom_line(aes(alpha = 0.1))
p1

p2 <- datos.gathered %>% filter(grupo == "Pu") %>%
  ggplot(aes(x=Wavelength,y=Intensity,color=archivo )) + 
  geom_line(aes(alpha = 0.1)) #parece ser que pu2 tiene problemas de calidad
p2

datos.gathered %>% filter(archivo == "Pu2") %>%
  ggplot(aes(x=Wavelength,y=Intensity,color= Nbarrido )) + 
  geom_line(aes(alpha = 0.1)) #el problema es específicamente con "Measure2"

datos.tidy <- datos.tidy%>% rownames_to_column('barrido') %>%
  filter(archivo != "Pu2" | Nbarrido != "Measure2") %>% #así que lo remuevo
  column_to_rownames('barrido')
datos.gathered <- filter(datos.gathered,barrido != "Pu2-B2") 

p3 <- datos.gathered %>% filter(grupo == "RD") %>%
  ggplot(aes(x=Wavelength,y=Intensity,color=archivo )) + 
  geom_line(aes(alpha = 0.1))
p3

p4 <- datos.gathered %>% filter(grupo == "TN") %>%
  ggplot(aes(x=Wavelength,y=Intensity,color=archivo )) + 
  geom_line(aes(alpha = 0.1))
p4

grid.arrange(p1,p2,p3,p4,ncol = 2, nrow = 2)#hay que graficar otra vez tras filtrar el espectro malo primero

#CONCLUSIón: se ve que los espectros no son tan diferentes como uno quisiera, probablemente generará algunos problemas a futuro
#¿quizás es necesario filtrar los espectros de otra forma?

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

#CONCLUSIÓN: los heatmaps de correlaciones no dan tanta información como se esperaba

#Análisis de componentes principales
df_pca <- datos.tidy %>% select(-grupo,-archivo,-Nbarrido) %>% prcomp(scale. = T) #la normalización no me permite escalar
df_out <- as.data.frame(df_pca$x)
labels <- datos.tidy %>% select(grupo,archivo,Nbarrido) %>% rownames_to_column('barrido') 
df_out <- df_out %>% rownames_to_column('barrido') %>% 
  left_join(labels, by = 'barrido') %>% column_to_rownames('barrido')

screeplot(df_pca)#las primeras dos componentes son suficientes

ggplot(df_out,aes(x=PC1,y=PC2,color=grupo )) + geom_point() #+ xlim()
#el PCA no parece generar un buen clustering, ¿será un problema con los datos o que la reducción dimensional no ayuda mucho?

#esto es para observar si hay algún efecto sobre el número de barrido, pero no parece ser el caso
ggplot(df_out,aes(x=PC1,y=PC2,color=Nbarrido )) + geom_point() #+ xlim()

#Pregunta: ¿puedo explicar alguna de las diferencias entre los grupos mediante alguna PC?
#si bien es cierto que el gráfico indica que probablemente no, alguna conclusión podría obtenerse
shapiro.test(df_out$PC1[df_out$grupo == "PE"])
shapiro.test(df_out$PC2[df_out$grupo == "PE"])
shapiro.test(df_out$PC1[df_out$grupo == "Pu"])
shapiro.test(df_out$PC2[df_out$grupo == "Pu"])
shapiro.test(df_out$PC1[df_out$grupo == "RD"])
shapiro.test(df_out$PC2[df_out$grupo == "RD"])
shapiro.test(df_out$PC1[df_out$grupo == "TN"])
shapiro.test(df_out$PC2[df_out$grupo == "TN"])#PC1 normal, PC2 no
bartlett.test(PC1~grupo,df_out)#varianzas homogéneas
aov_res <- aov(PC1 ~ grupo, data = df_out)
summary (aov_res)#PC1 no es suficientemente bueno para dividir entre grupos
TukeyHSD(aov_res)
kruskal.test(PC2 ~ grupo, data = df_out)#aunque PC2 parece ser interesante
pairwise.wilcox.test(df_out$PC2, df_out$grupo,p.adjust.method = "BH")

#CONCLUSIÓN: PC2 parece ser suficientemente decente para distinguir entre PE y resto, falla en separar otros pares

#otra pregunta: ¿hay alguna longitud de onda significativa? Podríamos analizar los loadings del PCA para determinarlo
#la idea es que si grafico loadings de PC1 vs. los de PC2 puedo observar longitudes de onda interesantes en aquellos lugares donde los puntos generados sean extremos
df_loadings <- as.data.frame(df_pca$rotation)
g <-df_loadings %>% 
  ggplot(aes(x=PC1, y =PC2)) +
  geom_point(alpha=0.4, size=1) + 
  geom_text(aes(label=rownames(df_loadings),hjust=0,vjust=0,
                alpha = ifelse(between(PC1,-0.026,0.11) & between(PC2,-0.09,0.18),0,1)))
g #¡hacer zoom para que sea más cómodo!
#resulta que las longitudes de onda que más agresivamente influyen en el PCA (y que puede que sean representativas) son 481.374,482.032 y 482.361 (corresponden a uno de los picos)
#puede que el filtrado sea muy agresivo o que convenga utilizar otros métodos para elegir longitudes de onda
#dicho eso, ahora puedo realizar ANOVAs o los equivalentes paramétricos
shapiro.test(datos.tidy$`481.374`[datos.tidy$grupo == "PE"])
shapiro.test(datos.tidy$`481.374`[datos.tidy$grupo == "Pu"])
shapiro.test(datos.tidy$`481.374`[datos.tidy$grupo == "RD"])
shapiro.test(datos.tidy$`481.374`[datos.tidy$grupo == "TN"])#este grupo no es normal
shapiro.test(datos.tidy$`482.032`[datos.tidy$grupo == "PE"])#no es normal
shapiro.test(datos.tidy$`482.032`[datos.tidy$grupo == "Pu"])
shapiro.test(datos.tidy$`482.032`[datos.tidy$grupo == "RD"])
shapiro.test(datos.tidy$`482.032`[datos.tidy$grupo == "TN"])#no es normal
shapiro.test(datos.tidy$`482.361`[datos.tidy$grupo == "PE"])
shapiro.test(datos.tidy$`482.361`[datos.tidy$grupo == "Pu"])
shapiro.test(datos.tidy$`482.361`[datos.tidy$grupo == "RD"])
shapiro.test(datos.tidy$`482.361`[datos.tidy$grupo == "TN"])#no es normal
kruskal.test(`481.374` ~ grupo, data = datos.tidy)#no es significativo
kruskal.test(`482.032` ~ grupo, data = datos.tidy)#no es significativo
kruskal.test(`482.361` ~ grupo, data = datos.tidy)#no es significativo
#pairwise.wilcox.test(datos.tidy$`481.374`, datos.tidy$grupo,p.adjust.method = "BH")

#extracción de características útiles para el análisis espectral
#recomendado por compañeros de trabajo
#lo utilizo para obtener el índice de una longitud de onda dada
featurevector <- datos.gathered %>% group_by(barrido) %>% 
  summarise(
            integral471a474.5 = trapz(Wavelength[between(Wavelength,471,474.5)],Intensity[between(Wavelength,471,474.5)]),
            integral488a497   = trapz(Wavelength[between(Wavelength,488,497  )],Intensity[between(Wavelength,488,497  )]),
            integral498a506   = trapz(Wavelength[between(Wavelength,498,506  )],Intensity[between(Wavelength,498,506  )]),
            integral509a515   = trapz(Wavelength[between(Wavelength,509,515  )],Intensity[between(Wavelength,509,515  )]),
            pico472.6  = interp1(Wavelength,Intensity,472.6 ),
            pico492.99 = interp1(Wavelength,Intensity,492.99),
            pico500.7  = interp1(Wavelength,Intensity,500.7 ),
            pico512.24 = interp1(Wavelength,Intensity,512.24)
            ) %>%
  column_to_rownames('barrido')

featurepca <- prcomp(featurevector,scale. = TRUE)
featureout <- as.data.frame(featurepca$x)
featureout <- featureout %>% rownames_to_column('barrido') %>% 
  left_join(labels, by = 'barrido') %>% column_to_rownames('barrido')

screeplot(featurepca)#las primeras dos componentes son suficientes

ggplot(featureout,aes(x=PC1,y=PC2,color=grupo )) + geom_point() #+ xlim(-1.5,1.5)
ggplot(featureout,aes(x=PC2,y=PC3,color=grupo )) + geom_point()
#parece ser un poco mejor el clustering pero no es suficiente; ¿quizás otras características funcionarían mejor?

#Analizo las PC
shapiro.test(featureout$PC1[featureout$grupo == "PE"])
shapiro.test(featureout$PC2[featureout$grupo == "PE"])
shapiro.test(featureout$PC1[featureout$grupo == "Pu"])
shapiro.test(featureout$PC2[featureout$grupo == "Pu"])
shapiro.test(featureout$PC1[featureout$grupo == "RD"])
shapiro.test(featureout$PC2[featureout$grupo == "RD"])
shapiro.test(featureout$PC1[featureout$grupo == "TN"])
shapiro.test(featureout$PC2[featureout$grupo == "TN"])#ambas normales
bartlett.test(PC1~grupo,featureout)
bartlett.test(PC2~grupo,featureout)#varianzas homogéneas
manova_res <- manova(cbind(PC1,PC2) ~ grupo, data = featureout)
summary(manova_res)
summary.aov(manova_res)#dividido por el número de test sólo es interesante PC2
aov_res <- aov(PC2 ~ grupo, data = featureout)
summary (aov_res)
TukeyHSD(aov_res)#las únicas diferencias significativas las veo para comparar TN con el resto

#CONCLUSIÓN: PC2 parece ser suficientemente decente para distinguir entre TN y resto, falla en separar otros pares

