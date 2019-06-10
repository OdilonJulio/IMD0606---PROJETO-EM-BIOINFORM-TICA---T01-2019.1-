
library(ggplot2)
library(tidyr)
library(dplyr)
library(lubridate)

setwd("~/IMD0606---PROJETO-EM-BIOINFORM-TICA---T01-2019.1-/IMD0606 - PROJETO EM BIOINFORMÁTICA - T01 (2019.1)")

# Expressão Gênica normalizada (Z-Score) do CD274 em tecidos Skin com câncer (Melanoma).
SkinExpGen <- read.csv('cBioPortal_data_SKIN_mRNA_Expression z-Scores(RNA Seq).csv')
LungExpGen <- read.csv('cBioPortal_data_LUNG_mRNA_Expression z-Scores(RNA Seq).csv')

# Ordenando pela Expressão Gênica em Z-Score
LungExpGen <- arrange(LungExpGen, desc(Exp_Z_CD274))
SkinExpGen <- arrange(SkinExpGen, desc(Exp_Z_CD274))

# Máximas Expressões Relativas
head(LungExpGen)
head(SkinExpGen)

max(LungExpGen$Exp_Z_CD274)
max(SkinExpGen$Exp_Z_CD274)

# Mínimas Expressões Relativas
tail(LungExpGen)
tail(SkinExpGen)

min(LungExpGen$Exp_Z_CD274)
min(SkinExpGen$Exp_Z_CD274)

# Histogramas e outros gráficos
hist(LungExpGen$Exp_Z_CD274, nclass=50)
hist(SkinExpGen$Exp_Z_CD274, nclass=50)

plot(dnorm(LungExpGen$Exp_Z_CD274))

hist(LungExpGen$Exp_Z_CD274, prob = T, nclass=50)
hist(SkinExpGen$Exp_Z_CD274, prob = T, nclass=50)

#x <- seq(min(LungExpGen$Exp_Z_CD274), max(LungExpGen$Exp_Z_CD274), 0.1)
plot(LungExpGen)

# Médias
mean(LungExpGen$Exp_Z_CD274)
mean(SkinExpGen$Exp_Z_CD274)

# Desvio Padrão
sd(LungExpGen$Exp_Z_CD274)
sd(SkinExpGen$Exp_Z_CD274)

#LungExpGenMax <- head(LungExpGen,100)
#LungExpGenMin <- tail(LungExpGen,100)

dim(SkinExpGen)
dim(LungExpGen)

## Filtros

# Calda dos muitos expressos
LungCaldaMax <- filter(LungExpGen, LungExpGen$Exp_Z_CD274 >= (sd(LungExpGen$Exp_Z_CD274)/2))
min(LungCaldaMax$Exp_Z_CD274)

# Calda dos poucos expressos
LungCaldaMin <- filter(LungExpGen, LungExpGen$Exp_Z_CD274 <= -(sd(LungExpGen$Exp_Z_CD274)/2))
max(LungCaldaMin$Exp_Z_CD274)

