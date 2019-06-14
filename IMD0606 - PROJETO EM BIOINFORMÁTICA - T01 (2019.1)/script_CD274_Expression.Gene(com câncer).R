



library(ggplot2)
library(swirl)
library(tidyr)
library(dplyr)
library(lubridate)
library(survival)

setwd("~/IMD0606---PROJETO-EM-BIOINFORM-TICA---T01-2019.1-/IMD0606 - PROJETO EM BIOINFORMÁTICA - T01 (2019.1)")

# Expressão Gênica normalizada (Z-Score) do CD274 em tecidos Skin com câncer (Melanoma).
SkinExpGen <- read.csv('cBioPortal_data_SKIN_mRNA_Expression z-Scores(RNA Seq).csv',  dec = ",")
LungExpGen <- read.csv('cBioPortal_data_LUNG_mRNA_Expression z-Scores(RNA Seq).csv',  dec = ",")

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

min(LungExpGen$Exp_Z_CD274)
min(SkinExpGen$Exp_Z_CD274)

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


hist(LungExpGen$Exp_Z_CD274, prob = T, nclass=50, main="Lung")
hist(SkinExpGen$Exp_Z_CD274, prob = T, nclass=50, main="Skin")

par(mfrow = c(1, 1))

### Filtros
# Selecionando pacientes ainda com câncer.
# Lung
LungExpGenCa <- filter(LungExpGen, LungExpGen$CANCER == "Recurred/Progressed")
dim(LungExpGenCa)
# Skin
SkinExpGenCa <- filter(SkinExpGen, SkinExpGen$CANCER == "Recurred/Progressed")
dim(SkinExpGenCa)


## Separando os indivíduos com maior e menor expressão.
# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, LungExpGenCa$Exp_Z_CD274 > median(LungExpGenCa$Exp_Z_CD274))
SkinExpGenCaMax <- filter(SkinExpGenCa, SkinExpGenCa$Exp_Z_CD274 > median(SkinExpGenCa$Exp_Z_CD274))
# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, LungExpGenCa$Exp_Z_CD274 < median(LungExpGenCa$Exp_Z_CD274))
SkinExpGenCaMin <- filter(SkinExpGenCa, SkinExpGenCa$Exp_Z_CD274 < median(SkinExpGenCa$Exp_Z_CD274))
  
# Superexpressos (principais pacientes)
summary(LungExpGenCaMax)
summary(SkinExpGenCaMax)

# Pouco expressos (principais pacientes)
summary(LungExpGenCaMin)
summary(SkinExpGenCaMin)

# Tratando a coluna SURVIVAL
LungExpGenCaMax <- mutate(LungExpGenCaMax, SURVIVAL=as.integer(SURVIVAL))
LungExpGenCaMin <- mutate(LungExpGenCaMin, SURVIVAL=as.integer(SURVIVAL))
SkinExpGenCaMax <- mutate(SkinExpGenCaMax, SURVIVAL=as.integer(SURVIVAL))
SkinExpGenCaMin <- mutate(SkinExpGenCaMin, SURVIVAL=as.integer(SURVIVAL))

LungExpGenCaMax <- arrange(LungExpGenCaMax, desc(SURVIVAL))
LungExpGenCaMin <- arrange(LungExpGenCaMin, desc(SURVIVAL))
SkinExpGenCaMax <- arrange(SkinExpGenCaMax, desc(SURVIVAL))
SkinExpGenCaMin <- arrange(SkinExpGenCaMin, desc(SURVIVAL))

## Survival Overall
SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

# SurvLung <- mutate(SurvLungMax, SURVIVALMin = SurvLungMin$SURVIVALMin); SurvLung
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge

fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)

autoplot(fit)
