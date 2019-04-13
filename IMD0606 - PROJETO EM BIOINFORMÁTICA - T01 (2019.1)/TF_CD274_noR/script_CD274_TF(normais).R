
library(ggplot2)
library(tidyr)
library(dplyr)
library(lubridate)



# Todos os TF do CD274 em tecidos normais com os respectivos R-Pearson em Lung e Skin
CD274_TF <- read.csv('PDL1.html.csv')

## Criando subset para Lung e Skin separadamente.

CD274_TF_Lung <- select(CD274_TF, tf_symbol, tf_score, r_lung, font)
head(CD274_TF_Lung)
CD274_TF_Skin <- select(CD274_TF, tf_symbol, tf_score, r_skin_not_sun_exposed_suprapubic, font)
head(CD274_TF_Skin)










summary(CD274_TF)
count(CD274_TF , tf_symbol)

