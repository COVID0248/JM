# ******************************************************************************************** #
# ******************************************************************************************** #
#                                              DATA                                            #
# ******************************************************************************************** #
# ******************************************************************************************** #
# Reading data in a short format
short_data <- read.csv("short_data.csv", sep = ";")

# Excluding NA
short_rmv.casos_sin_sintomas <- which(is.na(short_data$casos_sin_sintomas))
short_rmv.casos_activos <- which(is.na(short_data$casos_activos))
short_rmv.pcr <- which(is.na(short_data$pcr))
short_rmv.fallecidos <- which(is.na(short_data$fallecidos))

short_data$positividad <- (short_data$casos_con_sintomas + short_data$casos_sin_sintomas)/short_data$pcr

# Some transformations: Data for 100000 people
short_data$casos_con_sintomas <- (short_data$casos_con_sintomas/short_data$poblacion_region)*100000
short_data$casos_sin_sintomas <- (short_data$casos_sin_sintomas/short_data$poblacion_region)*100000
short_data$uci <- (short_data$uci/short_data$poblacion_region)*100000
short_data$pcr <- (short_data$pcr/short_data$poblacion_region)*100000

short_data$inmigrantes <- 100000*(short_data$inmigrantes/short_data$poblacion)
short_data$casos_activos <- 100000*(short_data$casos_activos/short_data$poblacion)
short_data$fallecidos <- 100000*(short_data$fallecidos/short_data$poblacion)
short_data$poblacion <- short_data$poblacion/100000

short_data$idse <- short_data$idse/1000 # IDSE from 0 to 1

# ************************************************************************ #
# Reading data in a long format
long_data <- read.csv("long_data.csv", sep = ";")

# Excluding NA
long_rmv.casos_sin_sintomas <- which(is.na(long_data$casos_sin_sintomas))
long_rmv.casos_activos <- which(is.na(long_data$casos_activos))
long_rmv.pcr <- which(is.na(long_data$pcr))
long_rmv.fallecidos <- which(is.na(long_data$fallecidos))

long_data$positividad <- (long_data$casos_con_sintomas + long_data$casos_sin_sintomas)/long_data$pcr

# Some transformations: Data for 100000 people
long_data$casos_con_sintomas <- (long_data$casos_con_sintomas/long_data$poblacion_region)*100000
long_data$casos_sin_sintomas <- (long_data$casos_sin_sintomas/long_data$poblacion_region)*100000
long_data$uci <- (long_data$uci/long_data$poblacion_region)*100000
long_data$pcr <- (long_data$pcr/long_data$poblacion_region)*100000

long_data$inmigrantes <- 100000*(long_data$inmigrantes/long_data$poblacion)
long_data$casos_activos <- 100000*(long_data$casos_activos/long_data$poblacion)
long_data$fallecidos <- 100000*(long_data$fallecidos/long_data$poblacion)
long_data$poblacion <- long_data$poblacion/100000

long_data$idse <- long_data$idse/1000 # IDSE from 0 to 1

# ******************************************************************************************** #
# ******************************************************************************************** #
#                                             MODEL                                            #
# ******************************************************************************************** #
# ******************************************************************************************** #
library(survival)
library(nlme)
library(JM)
library(joineRML)

max_tiempo <- max(long_data$tiempo2, short_data$tiempo2)
long_data$tiempo <- long_data$tiempo2/max_tiempo
short_data$tiempo <- short_data$tiempo2/max_tiempo

# Against the error in aeqSurv(Y): aeqSurv exception, an interval has effective length 0
short_data$tiempo3 <- short_data$tiempo + 0.001

long_data$log_casos_con_sintomas <- log(long_data$casos_con_sintomas+1e-07)
long_data$log_uci <- log(long_data$uci+1e-07)
long_data$log_casos_sin_sintomas <- log(long_data$casos_sin_sintomas+1e-07)
long_data$log_pcr <- log(long_data$pcr+1e-07)
long_data$log_casos_activos <- log(long_data$casos_activos+1e-07)
long_data$log_fallecidos <- log(long_data$fallecidos+1e-07)
long_data$log_positividad <- log(long_data$positividad+1e-07)

# ******************************************************************************************** #
# Univariate Joint Model 
# ******************************************************************************************** #
# uci

joint_fit_uci <- mjoint(formLongFixed = log_uci ~ tiempo,
                      formLongRandom = ~ tiempo | id,
                      formSurv = Surv(tiempo3, fim_cuarentena) ~ poblacion + densidad_poblacional + idse + 
                      hacinamiento + inmigrantes + capital_regional + capital_provincial + 
                      indice_ruralidad + aeropuerto + puerto + segunda_cuarentena,
                      data = long_data, survData = short_data, timeVar = "tiempo", verbose = F)
summary(joint_fit_uci)

# ******************************************************************************************** #
# casos_activos

long <- long_data[!is.na(long_data$casos_activos) & long_data$comuna != "Rengo" & 
                    long_data$comuna != "Punta Arenas",]
short <- short_data[!is.na(short_data$casos_activos),]

joint_fit_activos <- mjoint(formLongFixed = log_casos_activos ~ tiempo,
                      formLongRandom = ~ tiempo | id,
                      formSurv = Surv(tiempo3, fim_cuarentena) ~ poblacion + densidad_poblacional + idse + 
                      hacinamiento + inmigrantes + capital_regional + capital_provincial + 
                      indice_ruralidad + aeropuerto + puerto + segunda_cuarentena,
                      data = long, survData = short, timeVar = "tiempo", verbose = F)
summary(joint_fit_activos)

# ******************************************************************************************** #
# fallecidos - It is not significant

long <- long_data[!is.na(long_data$fallecidos) & long_data$comuna != "Angol" & 
                    long_data$comuna != "Antofagasta" & long_data$comuna != "Arica" & 
                    long_data$comuna != "Chillan" & long_data$comuna != "Chillan Viejo" &
                    long_data$comuna != "Hualpen" & long_data$comuna != "Las Condes" & 
                    long_data$comuna != "Lo Barnechea" & long_data$comuna != "Lonquimay" &
                    long_data$comuna != "Mejillones" & long_data$comuna != "Nueva Imperial" &
                    long_data$comuna != "Nunoa" & long_data$comuna != "Osorno" &
                    long_data$comuna != "Padre Las Casas" & long_data$comuna != "Providencia" & 
                    long_data$comuna != "Punta Arenas" & long_data$comuna != "Rengo" & 
                    long_data$comuna != "San Pedro de la Paz",]
short <- short_data[!is.na(short_data$fallecidos),]

joint_fit_fallecidos <- mjoint(formLongFixed = log_fallecidos ~ tiempo,
                      formLongRandom = ~ tiempo | id,
                      formSurv = Surv(tiempo3, fim_cuarentena) ~ poblacion + densidad_poblacional + idse + 
                        hacinamiento + inmigrantes + capital_regional + capital_provincial + 
                        indice_ruralidad + aeropuerto + puerto + segunda_cuarentena,
                      data = long, survData = short, timeVar = "tiempo", verbose = F)
summary(joint_fit_fallecidos)

# ******************************************************************************************** #
# positividad - It is not significant

long <- long_data[!is.na(long_data$positividad) & long_data$comuna != "Chillan" & 
                    long_data$comuna != "Chillan Viejo" & long_data$comuna != "Hualpen" & 
                    long_data$comuna != "Las Condes" & long_data$comuna != "Lo Barnechea" & 
                    long_data$comuna != "Nueva Imperial" & long_data$comuna != "Padre Las Casas" & 
                    long_data$comuna != "Providencia" & long_data$comuna != "San Pedro de la Paz" & 
                    long_data$comuna != "Vitacura",]
short <- short_data[!is.na(short_data$positividad),]

joint_fit_positivity <- mjoint(formLongFixed = log_positividad ~ tiempo,
                      formLongRandom = ~ tiempo | id,
                      formSurv = Surv(tiempo3, fim_cuarentena) ~ poblacion + densidad_poblacional + idse + 
                      hacinamiento + inmigrantes + capital_regional + capital_provincial + 
                      indice_ruralidad + aeropuerto + puerto + segunda_cuarentena,
                      data = long, survData = short, timeVar = "tiempo", verbose = F)
summary(joint_fit_positivity)

# ******************************************************************************************** #
# Bivariate Joint Model 
# ******************************************************************************************** #
# Coeficiente de Pearson
round(cor(long_data[,c("uci","casos_activos","positividad","fallecidos")], 
          method="pearson",use="na.or.complete"),2)

# ******************************************************************************************** #
# UCI + casos_activos - Both are significant
long <- long_data[!is.na(long_data$casos_activos) & long_data$comuna != "Rengo" & 
                    long_data$comuna != "Punta Arenas",]
short <- short_data[!is.na(short_data$casos_activos),]

joint_fit_uci_activos <- mjoint(formLongFixed = list("uci" = log_uci ~ tiempo,
                                             "casos_activos" = log_casos_activos ~ tiempo),
                        formLongRandom = list("uci" = ~ tiempo | id,
                                              "casos_activos" = ~ tiempo | id),
                        formSurv = Surv(tiempo3, fim_cuarentena) ~ poblacion + densidad_poblacional + idse + 
                        hacinamiento + inmigrantes + capital_regional + capital_provincial + 
                        indice_ruralidad + aeropuerto + puerto + segunda_cuarentena,
                        data = long, survData = short, timeVar = "tiempo", verbose = F)
summary(joint_fit_uci_activos) 

# ******************************************************************************************** #
# positividad + uci - uci is not significant
long <- long_data[!is.na(long_data$positividad) & long_data$comuna != "Chillan" & 
                    long_data$comuna != "Chillan Viejo" & long_data$comuna != "Hualpen" & 
                    long_data$comuna != "Las Condes" & long_data$comuna != "Lo Barnechea" & 
                    long_data$comuna != "Nueva Imperial" & long_data$comuna != "Padre Las Casas" & 
                    long_data$comuna != "Providencia" & long_data$comuna != "San Pedro de la Paz" & 
                    long_data$comuna != "Vitacura",]
short <- short_data[!is.na(short_data$casos_sin_sintomas),]

joint_fit_positividad_uci <- mjoint(formLongFixed = list("uci" = log_uci ~ tiempo,
                                             "positividad" = log_positividad ~ tiempo),
                        formLongRandom = list("uci" = ~ tiempo | id,
                                              "positividad" = ~ tiempo | id),
                        formSurv = Surv(tiempo3, fim_cuarentena) ~ poblacion + densidad_poblacional + idse + 
                          hacinamiento + inmigrantes + capital_regional + capital_provincial + 
                          indice_ruralidad + aeropuerto + puerto + segunda_cuarentena,
                        data = long, survData = short, timeVar = "tiempo", verbose = F)
summary(joint_fit_positividad_uci)

# ******************************************************************************************** #
# positividad + activos - Both are significant
long <- long_data[!is.na(long_data$positividad) & !is.na(long_data$casos_activos) &
                    long_data$comuna != "Rengo" & long_data$comuna != "Punta Arenas" &
                    long_data$comuna != "Chillan" & long_data$comuna != "Chillan Viejo" & 
                    long_data$comuna != "Hualpen" & long_data$comuna != "Las Condes" & 
                    long_data$comuna != "Lo Barnechea" & long_data$comuna != "Nueva Imperial" & 
                    long_data$comuna != "Padre Las Casas" & long_data$comuna != "Providencia" & 
                    long_data$comuna != "San Pedro de la Paz" & long_data$comuna != "Vitacura",]
short <- short_data[!is.na(short_data$casos_sin_sintomas) & !is.na(short_data$casos_activos) ,]

joint_fit_positividad_activos <- mjoint(formLongFixed = list("casos_activos" = log_casos_activos ~ tiempo,
                                             "positividad" = log_positividad ~ tiempo),
                        formLongRandom = list("casos_activos" = ~ tiempo | id,
                                              "positividad" = ~ tiempo | id),
                        formSurv = Surv(tiempo3, fim_cuarentena) ~ poblacion + densidad_poblacional + idse + 
                        hacinamiento + inmigrantes + capital_regional + capital_provincial + 
                        indice_ruralidad + aeropuerto + puerto + segunda_cuarentena,
                        data = long, survData = short, timeVar = "tiempo", verbose = F)
summary(joint_fit_positividad_activos) 

# ******************************************************************************************** #
# Trivariate Joint Model 
# ******************************************************************************************** #
# positividad + activos + UCI - UCI is not significant
long <- long_data[!is.na(long_data$positividad) & !is.na(long_data$casos_activos) &
                    long_data$comuna != "Rengo" & long_data$comuna != "Punta Arenas" &
                    long_data$comuna != "Chillan" & long_data$comuna != "Chillan Viejo" & 
                    long_data$comuna != "Hualpen" & long_data$comuna != "Las Condes" & 
                    long_data$comuna != "Lo Barnechea" & long_data$comuna != "Nueva Imperial" & 
                    long_data$comuna != "Padre Las Casas" & long_data$comuna != "Providencia" & 
                    long_data$comuna != "San Pedro de la Paz" & long_data$comuna != "Vitacura",]
short <- short_data[!is.na(short_data$casos_sin_sintomas) & !is.na(short_data$casos_activos) ,]

joint_fit_posi_activos_uci <- mjoint(formLongFixed = list("casos_activos" = log_casos_activos ~ tiempo,
                                             "uci" = log_uci ~ tiempo,
                                             "positividad" = log_positividad ~ tiempo),
                        formLongRandom = list("casos_activos" = ~ tiempo | id,
                                              "uci" = ~ tiempo | id,
                                              "positividad" = ~ tiempo | id),
                        formSurv = Surv(tiempo3, fim_cuarentena) ~ poblacion + densidad_poblacional + idse + 
                        hacinamiento + inmigrantes + capital_regional + capital_provincial + 
                        indice_ruralidad + aeropuerto + puerto + segunda_cuarentena,
                        data = long, survData = short, timeVar = "tiempo", verbose = F)
summary(joint_fit_posi_activos_uci)

# ******************************************************************************************** #
# Variable selection
# ******************************************************************************************** #
# positividad + activos
long <- long_data[!is.na(long_data$positividad) & !is.na(long_data$casos_activos) &
                    long_data$comuna != "Rengo" & long_data$comuna != "Punta Arenas" &
                    long_data$comuna != "Chillan" & long_data$comuna != "Chillan Viejo" & 
                    long_data$comuna != "Hualpen" & long_data$comuna != "Las Condes" & 
                    long_data$comuna != "Lo Barnechea" & long_data$comuna != "Nueva Imperial" & 
                    long_data$comuna != "Padre Las Casas" & long_data$comuna != "Providencia" & 
                    long_data$comuna != "San Pedro de la Paz" & long_data$comuna != "Vitacura",]
short <- short_data[!is.na(short_data$casos_sin_sintomas) & !is.na(short_data$casos_activos) ,]

# joint_fit_positividad_activos - densidad_poblacional (p-value = 0.83)
joint_fit_positividad_activos1 <- mjoint(formLongFixed = list("casos_activos" = log_casos_activos ~ tiempo,
                                             "positividad" = log_positividad ~ tiempo),
                        formLongRandom = list("casos_activos" = ~ tiempo | id,
                                              "positividad" = ~ tiempo | id),
                        formSurv = Surv(tiempo3, fim_cuarentena) ~ poblacion + idse + 
                        hacinamiento + inmigrantes + capital_regional + capital_provincial + 
                        indice_ruralidad + aeropuerto + puerto + segunda_cuarentena,
                        data = long, survData = short, timeVar = "tiempo", verbose = F)
summary(joint_fit_positividad_activos1) 

# joint_fit_positividad_activos1 - idse (p-value = 0.80)
joint_fit_positividad_activos2 <- mjoint(formLongFixed = list("casos_activos" = log_casos_activos ~ tiempo,
                                              "positividad" = log_positividad ~ tiempo),
                         formLongRandom = list("casos_activos" = ~ tiempo | id,
                                               "positividad" = ~ tiempo | id),
                         formSurv = Surv(tiempo3, fim_cuarentena) ~ poblacion +
                         hacinamiento + inmigrantes + capital_regional + capital_provincial + 
                         indice_ruralidad + aeropuerto + puerto + segunda_cuarentena,
                         data = long, survData = short, timeVar = "tiempo", verbose = F)
summary(joint_fit_positividad_activos2) 

# joint_fit_positividad_activos2 - aeropuerto (p-value = 0.73)
joint_fit_positividad_activos3 <- mjoint(formLongFixed = list("casos_activos" = log_casos_activos ~ tiempo,
                                              "positividad" = log_positividad ~ tiempo),
                         formLongRandom = list("casos_activos" = ~ tiempo | id,
                                               "positividad" = ~ tiempo | id),
                         formSurv = Surv(tiempo3, fim_cuarentena) ~ poblacion +
                         hacinamiento + inmigrantes + capital_regional + capital_provincial + 
                         indice_ruralidad + puerto + segunda_cuarentena,
                         data = long, survData = short, timeVar = "tiempo", verbose = F)
summary(joint_fit_positividad_activos3)

# joint_fit_positividad_activos3 - capital_provincial (p-value = 0.64)
joint_fit_positividad_activos4 <- mjoint(formLongFixed = list("casos_activos" = log_casos_activos ~ tiempo,
                                              "positividad" = log_positividad ~ tiempo),
                         formLongRandom = list("casos_activos" = ~ tiempo | id,
                                               "positividad" = ~ tiempo | id),
                         formSurv = Surv(tiempo3, fim_cuarentena) ~ poblacion + hacinamiento + 
                         inmigrantes + capital_regional + indice_ruralidad + puerto + 
                         segunda_cuarentena, data = long, survData = short, timeVar = "tiempo", 
                         verbose = F)
summary(joint_fit_positividad_activos4)

# joint_fit_positividad_activos4 - inmigrantes (p-value = 0.53)
joint_fit_positividad_activos5 <- mjoint(formLongFixed = list("casos_activos" = log_casos_activos ~ tiempo,
                                              "positividad" = log_positividad ~ tiempo),
                         formLongRandom = list("casos_activos" = ~ tiempo | id,
                                               "positividad" = ~ tiempo | id),
                         formSurv = Surv(tiempo3, fim_cuarentena) ~ poblacion + hacinamiento + 
                         capital_regional + indice_ruralidad + puerto + segunda_cuarentena, 
                         data = long, survData = short, timeVar = "tiempo", verbose = F)
summary(joint_fit_positividad_activos5)

# joint_fit_positividad_activos5 - capital_regional (p-value = 0.37)
joint_fit_positividad_activos6 <- mjoint(formLongFixed = list("casos_activos" = log_casos_activos ~ tiempo,
                                              "positividad" = log_positividad ~ tiempo),
                         formLongRandom = list("casos_activos" = ~ tiempo | id,
                                               "positividad" = ~ tiempo | id),
                         formSurv = Surv(tiempo3, fim_cuarentena) ~ poblacion + hacinamiento + 
                         indice_ruralidad + puerto + segunda_cuarentena, 
                         data = long, survData = short, timeVar = "tiempo", verbose = F)
summary(joint_fit_positividad_activos6)

# joint_fit_positividad_activos6 - puerto (p-value = 0.49)
joint_fit_positividad_activos7 <- mjoint(formLongFixed = list("casos_activos" = log_casos_activos ~ tiempo,
                                              "positividad" = log_positividad ~ tiempo),
                         formLongRandom = list("casos_activos" = ~ tiempo | id,
                                               "positividad" = ~ tiempo | id),
                         formSurv = Surv(tiempo3, fim_cuarentena) ~ poblacion + hacinamiento + 
                         indice_ruralidad + segunda_cuarentena, data = long, survData = short, 
                         timeVar = "tiempo", verbose = F)
summary(joint_fit_positividad_activos7)

# joint_fit_positividad_activos7 - segunda_cuarentena (p-value = 0.11)
joint_fit_positividad_activos8 <- mjoint(formLongFixed = list("casos_activos" = log_casos_activos ~ tiempo,
                                              "positividad" = log_positividad ~ tiempo),
                         formLongRandom = list("casos_activos" = ~ tiempo | id,
                                               "positividad" = ~ tiempo | id),
                         formSurv = Surv(tiempo3, fim_cuarentena) ~ poblacion + hacinamiento + 
                         indice_ruralidad, data = long, survData = short, timeVar = "tiempo", 
                         verbose = F)
summary(joint_fit_positividad_activos8)

# joint_fit_positividad_activos8 - poblacion (p-value = 0.12)
joint_fit_positividad_activos9 <- mjoint(formLongFixed = list("casos_activos" = log_casos_activos ~ tiempo,
                                              "positividad" = log_positividad ~ tiempo),
                         formLongRandom = list("casos_activos" = ~ tiempo | id,
                                               "positividad" = ~ tiempo | id),
                         formSurv = Surv(tiempo3, fim_cuarentena) ~ hacinamiento + indice_ruralidad, 
                         data = long, survData = short, timeVar = "tiempo", verbose = F)
summary(joint_fit_positividad_activos9); confint(joint_fit_positividad_activos9)
