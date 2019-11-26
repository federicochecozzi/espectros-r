read.csv.con.id <- function(.nombre.archivo, .carpeta){
  
  .archivo <- paste0(.carpeta, .nombre.archivo)  
  
  read_csv2(.archivo, skip = 1) %>% #la primera línea de los csv tiene texto redundante
    mutate(archivo = .nombre.archivo) %>% 
    mutate(archivo = str_remove(archivo, ".csv"))
}

cargar.muestras <- function(.nombre.archivos, .carpeta){
  
  archivos <- paste0(.carpeta, .nombre.archivos)
  
  map(.nombre.archivos, read.csv.con.id, .carpeta = .carpeta) %>% 
    bind_rows() %>% 
    mutate(archivo = factor(archivo))
  
}  