# Manipulación y flujo de trabajo general
library(tidyverse)

# Visualización adicional
library(grid)
library(gridExtra)
library(patchwork)

# Modelos estadísticos y pruebas
library(broom)
library(car)
library(lmtest)
library(MASS)
library(nnet)
library(pscl)
library(VGAM)
library(faraway)

# Evaluación de modelos
library(boot)
library(caret)
library(statmod)

# Lectura de archivos
library(readxl)

# Reportes y tablas
library(knitr)
library(kableExtra)

# Datos
malf <- read_xlsx("datos/malf.xlsx")

malf$presencia <- ifelse(!is.na(malf$talla_rn), 1, 0)
malf <- malf %>% 
  mutate(
    malf_num = as.numeric(malf_num),
    malf_mult = as.numeric(malf_mult),
    sexo_rn = as.factor(sexo_rn),
    comuna = as.factor(comuna),
    estacion = as.factor(estacion)
  )

malf <- malf %>%
  arrange(fecha_nac)

# Modelo Poisson
m_pois <- glm(malf_num ~ edad_madre + sexo_rn + comuna + estacion, 
              family = poisson, data = malf)
summary(m_pois)

# Modelo Quasi-Poisson
m_qpois <- glm(malf_num ~ edad_madre + sexo_rn + comuna + estacion, 
               family = quasipoisson, data = malf)
summary(m_qpois)

# Modelo Binomial
m_nb <- glm.nb(malf_num ~ edad_madre + sexo_rn + comuna + estacion, data = malf)
summary(m_nb)

# Modelo Zero Inflated Poisson
m_zip <- zeroinfl(malf_num ~ edad_madre + sexo_rn + comuna + estacion, 
                  data = malf, dist = "poisson")
summary(m_zip)

# Modelo Zero Inflated Negative Binomial
m_zinb <- zeroinfl(malf_num ~ edad_madre + sexo_rn + comuna + estacion, 
                   data = malf, dist = "negbin")
summary(m_zinb)


######

#
malf$w20_pm25_cs_10 <- malf$w20_pm25_cs/10
malf$w20_pm25_sp_10 <- malf$w20_pm25_sp/10

summary(malf$w20_pm25_cs)
summary(malf$w20_pm25_cs_10)

# Create new variable: pregnancy start date
malf$fecha_inicio <- malf$fecha_nac - (malf$edad_gest * 7)
malf$fecha_nac_re <-malf$fecha_nac

class(malf$fecha_nac)
class((malf$edad_gest))

# Load lubridate if not already loaded
library(lubridate)

# Calculate start of pregnancy
malf$fecha_inicio <- malf$fecha_nac - days(malf$edad_gest * 7)

# Optional: convert to Date class if you don't need time-of-day
malf$fecha_inicio <- as.Date(malf$fecha_inicio)

# Extract month
month_start <- month(malf$fecha_inicio)

# Create season variable
malf$season <- case_when(
  month_start %in% c(12, 1, 2) ~ "summer",
  month_start %in% c(3, 4, 5)  ~ "fall",
  month_start %in% c(6, 7, 8)  ~ "winter",
  month_start %in% c(9, 10, 11) ~ "spring",
  TRUE ~ NA_character_
)

malf$season_cold <- ifelse(malf$season %in% c("fall", "winter"), 1, 0)
malf$season_winter_only <- ifelse(malf$season == "winter", 1, 0)


###Among all
#Using central site
malf_log<- glm(presencia~w20_pm25_cs_10 
               + edad_madre + sexo_rn + a_nac + season
               , family="binomial" (link="logit"), data=malf)
summary(malf_log)
exp(malf_log$coefficients[-1])
exp(confint.default(malf_log)) 

#using spatial
malf_log_sp<- glm(presencia~w20_pm25_sp_10
                  + edad_madre + sexo_rn + a_nac + season
                  , family="binomial" (link="logit"), data=malf)
summary(malf_log_sp)
exp(malf_log_sp$coefficients[-1])
exp(confint.default(malf_log_sp)) 
####No relationship between exposure in first 20 weeks and risk of malformation

###SEEING IF SEASON OF CONCEPTION MAKES A DIFFERENCE

malf_winter_fall<- glm(presencia~season_cold 
                  + edad_madre + sexo_rn + a_nac
                  , family="binomial" (link="logit"), data=malf)
summary(malf_winter_fall)
exp(malf_winter_fall$coefficients[-1])
exp(confint.default(malf_winter_fall)) 
#No significant increased risk of malformation with conception in winter or fall

malf_winter<- glm(presencia~season_winter_only 
                       + edad_madre + sexo_rn + a_nac
                       , family="binomial" (link="logit"), data=malf)
summary(malf_winter)
exp(malf_winter$coefficients[-1])
exp(confint.default(malf_winter)) 
#No significant increased risk of malformation with conception in winter


###SEPARATE MODEL FOR ONLY THOSE CONCEIVED IN WINTER 
winter_data <- subset(malf, season_winter_only == 1)

####
library(gmodels)

table(winter_data$presencia, useNA = "always")

summary(winter_data$t1_pm25_cs[which(winter_data$presencia==0)])
summary(winter_data$t1_pm25_cs[which(winter_data$presencia==1)])

summary(winter_data$t1_pm25_sp[which(winter_data$presencia==0)])
summary(winter_data$t1_pm25_sp[which(winter_data$presencia==1)])

summary(winter_data$w20_pm25_cs[which(winter_data$presencia==0)])
summary(winter_data$w20_pm25_cs[which(winter_data$presencia==1)])



####Among those conceived in Winter
#Using central site
malf_log<- glm(presencia~w20_pm25_cs 
             + edad_madre + sexo_rn + a_nac
             , family="binomial" (link="logit"), data=winter_data)
summary(malf_log)
exp(malf_log$coefficients[-1])
exp(confint.default(malf_log)) 

#using spatial
malf_log_sp<- glm(presencia~w20_pm25_sp
                  + edad_madre + sexo_rn + a_nac
                  , family="binomial" (link="logit"), data=winter_data)
summary(malf_log_sp)
exp(malf_log_sp$coefficients[-1])
exp(confint.default(malf_log_sp)) 
###no relationship between exposure in the first 20 weeks and risk of malformation

####
summary(malf$w20_pm25_cs[which(malf$season_winter_only==0)])
summary(malf$w20_pm25_cs[which(malf$season_winter_only==1)])
#Clear indications of higher levels of PM according to the central site 
#among those who were conceived in winter versus not

summary(malf$w20_pm25_sp[which(malf$season_winter_only==0)])
summary(malf$w20_pm25_sp[which(malf$season_winter_only==1)])
#Clear indications of higher levels of PM according to the LUR
#among those who were conceived in winter versus not


###descriptive data for resumen
table(malf$presencia, useNA = "always")
summary(malf$w20_pm25_cs)
summary(malf$w20_pm25_sp)
