
#setwd("D:/UserData/BaseARG/Perfiles")

setwd("~/Documents/R")
perfiles <- read.csv("perfiles_2014-07-01.csv")#[,c(-1,-2)])

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

# How many calicata have same number/name at each data frame?
table(as.data.frame(table(locations$NAME))$Freq)
table(as.data.frame(table(calicatas$numero))$Freq)

#merge data frames
coord.profile <- merge(x = locations, y = calicatas, by.x = "NAME", by.y = "numero", all.x=TRUE, all.y=T)


# comparison between mosaico from coordinates and mosaico from sisinta 
coord.profile$mosaic_match <- as.numeric(coord.profile$mosaico_si == coord.profile$mosaico_qg)
coord.profile$id <- gsub(pattern='^', replacement="http://sisinta.inta.gov.ar/perfiles/" , x=coord.profile$id, fixed= F)  
# csv of matches and no-matches 
write.csv(coord.profile,"coord.profile.jul14.csv")


#### join coordinates to profiles
setwd("~/Documents/R")
coord <- read.csv("coord.profile.jul14.389p.csv")
data <- read.csv("perfiles_2014-07-01.csv")

calib.data <- merge(x = coord,y = data, by.x = "id", by.y = "perfil_id", all.x = T)

length(unique(calib.data$id))
