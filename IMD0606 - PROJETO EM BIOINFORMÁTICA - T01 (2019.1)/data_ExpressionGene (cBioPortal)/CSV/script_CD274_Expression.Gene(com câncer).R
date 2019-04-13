
library(ggplot2)
library(tidyr)
library(dplyr)
library(lubridate)



# Expressão Gênica do CD274 em tecidos Skin com câncer (Melanoma).
SkinExpGen <- read.csv('cBioPortal_data_Skin.csv')
LungExpGen <- read.csv('cBioPortal_data_Lung.csv')

# Ordenando pela Expressão Gênica
LungExpGen <- arrange(LungExpGen, desc(CD274))
SkinExpGen <- arrange(SkinExpGen, desc(CD274))

head(LungExpGen)
head(SkinExpGen)

dim(SkinExpGen)
dim(LungExpGen)

help(arrange)







summary(CD274_TF)
count(CD274_TF , tf_symbol)

