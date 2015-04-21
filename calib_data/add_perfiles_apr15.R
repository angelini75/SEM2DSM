
#setwd("D:/UserData/BaseARG/Perfiles")

setwd("~/Documents/R")
perfiles <- read.csv("perfiles_2015-04-17.csv")#[,c(-1,-2)])

# # inserting a column between "perfil_observaciones" and "perfil_ubicacion_descripcion"
# perfiles$perfil_capacidad_clase[perfiles$perfil_capacidad_clase != ""]
# a <- perfiles[perfiles$perfil_capacidad_clase != "",]
# write.csv(a, "a.csv")
# #### ATENTION!!!
# # OPENING A.CSV WITH LIBREOFICE (CALC)  and drag data leaving "perfil_capacidad_clase" field empty 
# a <- read.csv("a.csv")#
# rownames(a)<-a$X
# a <- a[,-1]
# perfiles <- perfiles[perfiles$perfil_capacidad_clase == "",]
# perfiles <- rbind(perfiles,a)

# from perfiles (with horizons) to calicatas (without orizons)
#names(perfiles)

#as.data.frame(table(perfiles$perfil_ubicacion_mosaico))

calicatas <- unique(data.frame("id"=perfiles$perfil_id, "numero"=perfiles$perfil_numero, 
                            "xy"=perfiles$perfil_ubicacion_coordenadas, "fecha"=perfiles$perfil_fecha,
                            "mosaico_si"=perfiles$perfil_ubicacion_mosaico))
calicatas$fecha <- as.POSIXlt.date(calicatas$fecha)



# delete profiles newer than year 2000
calicatas <- calicatas[calicatas$fecha < 2000-01-01,]

# delete profiles without perfil_numero (numero) & 
calicatas <- calicatas[calicatas$numero != "",]
#calicatas[345,]c(5,3)] <- "" # columns with errors
#calicatas_noxy <-calicatas[calicatas$xy == "",]






# standarise profile names
calicatas$numero <- gsub(" C", "" , calicatas$numero, fixed=TRUE)
calicatas$numero <- gsub(" (FAO-UNESCO)", "" , calicatas$numero, fixed=TRUE)
calicatas$numero <- gsub(" FAO", "" , calicatas$numero, fixed=TRUE)
calicatas$numero <- gsub(" (FAO)", "" , calicatas$numero, fixed=TRUE)
calicatas$numero <- gsub("-C", "" , calicatas$numero, fixed=TRUE)
calicatas$numero <- gsub(" ", "" , calicatas$numero, fixed=TRUE)
calicatas$numero <- gsub("C", "" , calicatas$numero, fixed=TRUE)
calicatas$numero <- gsub("/", "-" , calicatas$numero, fixed=TRUE)
calicatas$numero <- gsub("RP", "" , calicatas$numero, fixed=TRUE)

#write.csv(calicatas,"calicatas.csv")

# load locations (x,y): two sources: cali cruzate and santiago points (called 123)
locations <- read.csv("locations_20_apr_hojas.csv")
locations <-locations[locations$NAME != "",]
names(locations)[5] <- "mosaico_qg"

# standarise profile names
locations$NAME <- gsub(" C", "" , locations$NAME, fixed=TRUE)
locations$NAME <- gsub("-C", "" , locations$NAME, fixed=TRUE)
locations$NAME <- gsub("C-", "" , locations$NAME, fixed=TRUE)
locations$NAME <- gsub(" C", "" , locations$NAME, fixed=TRUE)
locations$NAME <- gsub("C", "" , locations$NAME, fixed=TRUE)
locations$NAME <- gsub(" ", "" , locations$NAME, fixed=TRUE)
locations$NAME <- gsub("/", "-" , locations$NAME, fixed=TRUE)


# Remove "0" as first character
locations$NAME <- gsub(pattern='^0', replacement="" , x=locations$NAME,fixed= F)
locations$NAME <- gsub(pattern='^-0', replacement="" , x=locations$NAME,fixed= F)
calicatas$numero <- gsub(pattern='^0', replacement="" , x=calicatas$numero, fixed= F)

#Remove: 0 from mosaico names ????-0?-0? | words from mosaicos sisinta | white spaces
locations$mosaico_qg <- gsub(pattern='-0', replacement="-" , x=locations$mosaico_qg, fixed= F)
calicatas$mosaico_si <- gsub(pattern='-0', replacement="-" , x=calicatas$mosaico_si, fixed= F)
calicatas$mosaico_si <- gsub(pattern=' -', replacement="-" , x=calicatas$mosaico_si, fixed= F)
calicatas$mosaico_si <- gsub(pattern='- ', replacement="-" , x=calicatas$mosaico_si, fixed= F)  
calicatas$mosaico_si <- gsub(pattern=' \\w+', replacement="" , x=calicatas$mosaico_si, fixed= F)  
calicatas$mosaico_si <- gsub(pattern=' ', replacement="" , x=calicatas$mosaico_si, fixed= F)  
calicatas$mosaico_si
calicatas$mosaico_si[488]<-"3360-26-3"
# Selection profiles in mosaico 3560-15:17-1:4
ca <- calicatas[NULL,]
for(i in 14:17){
  for(j in 1:4){
   ca<- rbind(ca,calicatas[calicatas$mosaico_si==paste("3560",i,j, sep="-"),])
  }
}
ca<-ca[ca$xy!="",]

# How many calicata have same number/name at each data frame?
table(as.data.frame(table(locations$NAME))$Freq)
table(as.data.frame(table(calicatas$numero))$Freq)
n<-as.data.frame(table(calicatas$numero))
n[n$Freq>1,]

#merge data frames
coord.profile <- merge(x = locations, y = calicatas, by.x = "NAME", by.y = "numero", all.x=TRUE, all.y=T)


# comparison between mosaico from coordinates and mosaico from sisinta 
coord.profile$mosaic_match <- as.numeric(coord.profile$mosaico_si == coord.profile$mosaico_qg)
coord.profile$id <- gsub(pattern='^', replacement="http://sisinta.inta.gov.ar/perfiles/" , x=coord.profile$id, fixed= F)  
# csv of matches and no-matches 
write.csv(coord.profile,"coord.profile.apr15.csv")


#merge coord.profile with bea (locations without profiles)
bea <- read.csv("bea.csv")

m <- merge(bea, coord.profile,by="NAME", all.x=T, all.y=F)
m <- m[!is.na(m$id.y),]
m <- m[,c(14,4)]
names(m)<- c("id", "xy")
#### join coordinates to profiles
# setwd("~/Documents/R")
#coord <- read.csv("coord.profile.jul14.389p.csv")
coord <- ca[,c(1,3)]
data <- read.csv("perfiles_2015-04-17.csv")

calib.data <- merge(x = coord,y = data, by.x = "id", by.y = "perfil_id", all.x = T)
calib.data<- calib.data[!is.na(calib.data$analitico_registro),]
length(unique(calib.data$id))

### agregar script para adosar datos de calib.data a calib.data-0.1.1.csv

header <- c("a_registro","a_humedad","a_s","a_t","a_ph_pasta","a_ph_h2o","a_ph_kcl",
            "a_resistencia_pasta","a_base_ca","a_base_mg","a_base_k","a_base_na",
            "a_arcilla","a_materia_organica_c","a_materia_organica_n","a_limo_2_20",
            "a_limo_2_50","a_arena_muy_fina","a_arena_fina","a_arena_media","a_arena_gruesa",
            "a_arena_muy_gruesa","a_ca_co3","a_agua_15_atm","a_agua_util","a_conductividad",
            "a_h","a_saturacion_t","a_saturacion_s_h","a_densidad_aparente",
            "a_profundidad_muestra","a_agua_3_atm","a_materia_organica_cn",
            "barnices","co3","color_humedo_hvc","color_seco_hvc","concreciones",
            "consistencia_en_seco","consistencia_en_humedo","consistencia_adhesividad",
            "consistencia_plasticidad","estructura_tipo","estructura_clase",
            "estructura_grado","formaciones_especiales","humedad","id","limite_tipo",
            "limite_forma","moteados","p_id","p_profundidad_napa",           #"p_numero",
            "p_cobertura_vegetal","p_material_original","p_modal","p_fecha",
            "p_vegetacion_o_cultivos","p_observaciones","p_capacidad_clase",
            "p_u_descripcion","p_u_recorrido","p_u_mosaico",
            "p_u_aerofoto","p_paisaje_tipo",                                 #"p_u_coordenadas",
            "p_paisaje_forma","p_paisaje_simbolo","p_humedad_clase","p_erosion_clase",
            "p_erosion_subclase","p_pedregosidad_clase","p_pedregosidad_subclase",
            "p_serie_nombre","p_serie_simbolo","p_serie_descripcion","p_grupo_codigo",
            "p_grupo_descripcion","p_fase_codigo","p_fase_nombre","p_sales",
            "p_uso_de_la_tierra","p_relieve","p_permeabilidad","p_escurrimiento",
            "p_pendiente",,"p_anegamiento","p_drenaje","p_posicion","ph",
            "profundidad_inferior","profundidad_superior","raices","textura","tipo")

options(digits = 8) 
q <-as.double(unlist(strsplit(gsub(")","", gsub("POINT [(]", "", as.character(calib.data$xy))), split = " ")))
calib.data$X <- q[seq(1,length(q), 2)]
calib.data$Y <- q[seq(0,length(q), 2)]


names(calib.data) <- header
head(calib.data)
setwd("/media/L0135974_DATA/UserData/BaseARG/2_Calibration")
D <- read.csv("calib.data-0.2.csv")
length(unique(D$X))
n <- as.data.frame(table(D$p_id))
nodata <- as.character(n[n$Freq==1,][,1])
#any(D$p_id != c("11/834 C", "12/778-C"))

D[D$p_id== c("06/118 C","11/834 C","12/778-C"),]
D1 <- D[D$p_id!= c("12/778-C","11/834 C"),]
