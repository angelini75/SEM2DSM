#### This script clean up the calibration data
# importing calibration data


# header <- c("a_registro","a_humedad","a_s","a_t","a_ph_pasta","a_ph_h2o","a_ph_kcl",
#             "a_resistencia_pasta","a_base_ca","a_base_mg","a_base_k","a_base_na",
#             "a_arcilla","a_materia_organica_c","a_materia_organica_n","a_limo_2_20",
#             "a_limo_2_50","a_arena_muy_fina","a_arena_fina","a_arena_media","a_arena_gruesa",
#             "a_arena_muy_gruesa","a_ca_co3","a_agua_15_atm","a_agua_util","a_conductividad",
#             "a_h","a_saturacion_t","a_saturacion_s_h","a_densidad_aparente",
#             "a_profundidad_muestra","a_agua_3_atm","a_materia_organica_cn",
#             "barnices","co3","color_humedo_hvc","color_seco_hvc","concreciones",
#             "consistencia_en_seco","consistencia_en_humedo","consistencia_adhesividad",
#             "consistencia_plasticidad","estructura_tipo","estructura_clase",
#             "estructura_grado","formaciones_especiales","humedad","id","limite_tipo",
#             "limite_forma","moteados","p_id","p_profundidad_napa",           #"p_numero",
#             "p_cobertura_vegetal","p_material_original","p_modal","p_fecha",
#             "p_vegetacion_o_cultivos","p_observaciones","p_capacidad_clase",
#             "p_u_descripcion","p_u_recorrido","p_u_mosaico",
#             "p_u_aerofoto","p_paisaje_tipo",                                 #"p_u_coordenadas",
#             "p_paisaje_forma","p_paisaje_simbolo","p_humedad_clase","p_erosion_clase",
#             "p_erosion_subclase","p_pedregosidad_clase","p_pedregosidad_subclase",
#             "p_serie_nombre","p_serie_simbolo","p_serie_descripcion","p_grupo_codigo",
#             "p_grupo_descripcion","p_fase_codigo","p_fase_nombre","p_sales",
#             "p_uso_de_la_tierra","p_relieve","p_permeabilidad","p_escurrimiento",
#             "p_pendiente","p_anegamiento","p_drenaje","p_posicion","ph",
#             "profundidad_inferior","profundidad_superior","raices","textura","tipo")
# names(D) <- header
#### Working directory
setwd("/media/L0135974_DATA/UserData/BaseARG/2_Calibration")
####Load data
D <- read.csv("calib.data-0.2.csv")

#Counting number of profiles
length(unique(D$X))
#delete profiles with 1 horizon
n <- as.data.frame(table(D$p_id))
nodata <- as.character(n[n$Freq==1,][,1])
#any(D$p_id != c("11/834 C", "12/778-C"))
####### it fails
D[D$p_id== c("06/118 C","11/834 C","12/778-C"),]
D1 <- D[D$p_id!= c("12/778-C","11/834 C"),]
