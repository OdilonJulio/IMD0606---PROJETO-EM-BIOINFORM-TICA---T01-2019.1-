
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

dim(SkinExpGen)
dim(LungExpGen)

max(LungExpGen$Exp_Z_CD274)
max(SkinExpGen$Exp_Z_CD274)

# Mínimas Expressões Relativas
tail(LungExpGen)
tail(SkinExpGen)

min(LungExpGen$Exp_Z_CD274)
min(SkinExpGen$Exp_Z_CD274)

# Médidas de dispersão e de tendência central.
sd(LungExpGen$Exp_Z_CD274)
sd(SkinExpGen$Exp_Z_CD274)

mean(LungExpGen$Exp_Z_CD274)
mean(SkinExpGen$Exp_Z_CD274)

var(LungExpGen$Exp_Z_CD274)
var(SkinExpGen$Exp_Z_CD274)

median(LungExpGen$Exp_Z_CD274)
median(SkinExpGen$Exp_Z_CD274)

# Histogramas e outros gráficos
par(mfrow = c(1, 2)) # Possibilita visualizar vários plots ao mesmo tempo.
hist(LungExpGen$Exp_Z_CD274, nclass=500, main="Lung")
hist(SkinExpGen$Exp_Z_CD274, nclass=500, main="Skin")

qqnorm(LungExpGen$Exp_Z_CD274, main="Lung")
qqline(LungExpGen$Exp_Z_CD274, col='red')
qqnorm(SkinExpGen$Exp_Z_CD274, main="Skin")
qqline(SkinExpGen$Exp_Z_CD274, col='red')

plot(LungExpGen, main="Lung")
plot(SkinExpGen, main="Skin")

hist(LungExpGen$Exp_Z_CD274, prob = T, nclass=50, main="Lung")
hist(SkinExpGen$Exp_Z_CD274, prob = T, nclass=50, main="Skin")

par(mfrow = c(1, 1))

## Filtros
## Separando os indivíduos com maior e menor expressão em 
## 2,5% em cada uma das extremidades.

# Superexpressos.
LungExpGenMax <- head(LungExpGen,0.025*length(LungExpGen$Exp_Z_CD274))
SkinExpGenMax <- head(SkinExpGen,0.025*length(SkinExpGen$Exp_Z_CD274))

# Pouco expressos.
LungExpGenMin <- tail(LungExpGen,0.025*length(LungExpGen$Exp_Z_CD274))
SkinExpGenMin <- tail(SkinExpGen,0.025*length(SkinExpGen$Exp_Z_CD274))

# Avaliando distâncias.
# Distância entre o menor dos maiores e o maior dos menores - Lung
abs(min(LungExpGenMax$Exp_Z_CD274)/max(LungExpGenMin$Exp_Z_CD274))

# Distância entre o menor dos maiores e o maior dos menores - Skin
abs(min(SkinExpGenMax$Exp_Z_CD274)/max(SkinExpGenMin$Exp_Z_CD274))



## COLETAR A SOBREVIDA DESSES INDIVÍDUOS