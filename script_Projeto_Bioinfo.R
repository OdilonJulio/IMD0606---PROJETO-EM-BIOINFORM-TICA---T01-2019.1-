
########################################################################################
############################# Disciplina: Projeto em Bioinformática ####################
############################# Professores: Jorge e César ###############################
############################# Orientadora: Professora Tirzah ###########################
############################# Aluno: Odilon Júlio dos Santos ###########################
########################################################################################

# Bibliotecas necessárias
library(ggplot2)
library(swirl)
library(tidyr)
library(dplyr)
library(lubridate)
library(survival)
library(stringr)

### Acessando e analisando o banco de dados com os TFs relacionados ao CD274 nos tecidos Lung e Skin.

setwd("~/IMD0606---PROJETO-EM-BIOINFORM-TICA---T01-2019.1-")

TFs <- read_csv("PDL1.html.csv")

View(TFs)

dim(TFs) # Existem 825 TFs relacionados ao CD274

class(TFs)

head(TFs)

names(TFs)

#################################################################################################
### O objetivo agora é selecionar aqueles fatores que transcrição que mais fortemente influenciam
### na expressão do gene CD274.

## LUNG
# Selecionando os fatores de transcrição que possuem R-Pearson maior que 0.6 e menor que -0.4.

TFsRLungMax <- filter(TFs, r_lung > 0.5)
dim(TFsRLungMax)
head(TFsRLungMax)
#View(TFsRLungMax)

TFsRLungMin <- filter(TFs, r_lung < -0.4)
dim(TFsRLungMin)
head(TFsRLungMin)
#View(TFsRLungMin)

## SKIN
# Selecionando os fatores de transcrição que possuem R-Pearson maior que 0.5 e menor que -0.3.

TFsRSkinMax <- filter(TFs, r_skin_sun_exposed_lower_leg > 0.5)
dim(TFsRSkinMax)
head(TFsRSkinMax)
#View(TFsRSkinMax)

TFsRSkinMin <- filter(TFs, r_skin_sun_exposed_lower_leg < -0.3)
dim(TFsRSkinMin)
head(TFsRSkinMin)
#View(TFsRSkinMin)

#################################################################################################
### Agora, iremos montar datasets distintos para Lung e para Skin

## SKIN
# Construindo datasets com o nome do TF e sua correlação com o CD274
TFsRLungMax <- select(TFsRLungMax, tf_symbol, r_lung)
dim(TFsRLungMax)
head(TFsRSkinMax)

TFsRLungMax <- unique(TFsRLungMax) # Retirando duplicatas. 
dim(TFsRLungMax)
head(TFsRLungMax)

TFsRLungMin <- select(TFsRLungMin, tf_symbol, r_lung)
dim(TFsRLungMin)
head(TFsRLungMin)

TFsRLungMin <- unique(TFsRLungMin) # Retirando duplicatas. 
dim(TFsRLungMin)
head(TFsRLungMin)

## SKIN
# Construindo datasets com o nome do TF e sua correlação com o CD274
TFsRSkinMax <- select(TFsRSkinMax, tf_symbol, r_skin_sun_exposed_lower_leg)
dim(TFsRSkinMax)
head(TFsRSkinMax)

TFsRSkinMax <- unique(TFsRSkinMax) # Retirando duplicatas. 
dim(TFsRSkinMax)
head(TFsRSkinMax)

TFsRSkinMin <- select(TFsRSkinMin, tf_symbol, r_skin_sun_exposed_lower_leg)
dim(TFsRSkinMin)
head(TFsRSkinMin)

TFsRSkinMin <- unique(TFsRSkinMin) # Retirando duplicatas. 
dim(TFsRSkinMin)
head(TFsRSkinMin)

#################################################################################################
### Quais foram os TFs econtrados? 

## LUNG

# Correlação positiva
TFsRLungMax$tf_symbol
#View(TFsRLungMax$tf_symbol)
# Resultado:
#[1] "IRF1"    "RELA"    "TEAD4"   "CREM"    "ETS1"    "ETS2"    "BACH1"   "HMGA1"   "HSF1"    "NFKB2"  
#[11] "NKX3-1"  "RARA"    "REL"     "RELB"    "STAT5A"  "HIVEP2"  "PML"     "BCL3"    "NFYA"    "GATAD2A"
#[21] "TEAD3"

# Correlação negativa
TFsRLungMin$tf_symbol
#View(TFsRLungMin$tf_symbol)
# Resultado:
#[1] "PATZ1"  "ZNF554" "THAP11"


## SKIN

# Correlação positiva
TFsRSkinMax$tf_symbol
#View(TFsRSkinMax$tf_symbol)
# Resultado:
#[1] "STAT1" "IKZF1" "STAT2" "BATF3" "IRF7"  "PML"

# Correlação negativa
TFsRSkinMin$tf_symbol
#View(TFsRSkinMin$tf_symbol)
# Resultado:
#[1] "RARG" "RXRA" "TEF" 

#################################################################################################
### Já sabemos quais são os fatores de transcrição procurados. Devemos, agora, encontrar a 
### expressão gênica deles, preferencialmente, nos mesmos pacientes. Tal dado, foi retirado do
### banco TCGA, por meio do cBioPortal.org

# Expressão Gênica normalizada (Z-Score) do CD274 em tecidos Skin com câncer (Melanoma) e em 
# tecidos Lung com câncer (Adenocarcinoma).

LungExpGen <- read_csv("cBioPortal_data_LUNG_mRNA_Expression z-Scores(RNA Seq).csv")
SkinExpGen <- read_csv("cBioPortal_data_SKIN_mRNA_Expression z-Scores(RNA Seq).csv")
View(LungExpGen)
View(SkinExpGen)

# Ordenando pela Expressão Gênica em Z-Score
LungExpGen <- arrange(LungExpGen, desc(Exp_Z_CD274))
SkinExpGen <- arrange(SkinExpGen, desc(Exp_Z_CD274))

# Máximas Expressões Relativas
head(LungExpGen$Exp_Z_CD274)
head(SkinExpGen$Exp_Z_CD274)

dim(SkinExpGen)
dim(LungExpGen)

max(LungExpGen$Exp_Z_CD274)
max(SkinExpGen$Exp_Z_CD274)

min(LungExpGen$Exp_Z_CD274)
min(SkinExpGen$Exp_Z_CD274)

# Mínimas Expressões Relativas
tail(LungExpGen$Exp_Z_CD274)
tail(SkinExpGen$Exp_Z_CD274)

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

#jpeg("plot_histogramas.jpeg", pointsize = 20, width = 800, height = 800)
par(mfrow = c(1, 2)) # Possibilita visualizar vários plots ao mesmo tempo.
hist(LungExpGen$Exp_Z_CD274, nclass=500, main="Lung")
hist(SkinExpGen$Exp_Z_CD274, nclass=500, main="Skin")
#dev.off()

#jpeg("plot_testeNorm.jpeg", pointsize = 20, width = 800, height = 800)
par(mfrow = c(1, 2))
qqnorm(LungExpGen$Exp_Z_CD274, main="Lung")
qqline(LungExpGen$Exp_Z_CD274, col='red')
qqnorm(SkinExpGen$Exp_Z_CD274, main="Skin")
qqline(SkinExpGen$Exp_Z_CD274, col='red')
#dev.off()

#jpeg("plot_histograma2.jpeg", pointsize = 20, width = 800, height = 800)
par(mfrow = c(1, 2))
hist(LungExpGen$Exp_Z_CD274, prob = T, nclass=50, main="Lung")
hist(SkinExpGen$Exp_Z_CD274, prob = T, nclass=50, main="Skin")
#dev.off()

class(LungExpGen$BACH1)

par(mfrow = c(1, 1))

### Filtros
# Selecionando pacientes ainda com câncer.
# Lung
LungExpGenCa <- filter(LungExpGen, LungExpGen$CANCER == "Recurred/Progressed")
dim(LungExpGenCa)
# Skin
SkinExpGenCa <- filter(SkinExpGen, SkinExpGen$CANCER == "Recurred/Progressed")
dim(SkinExpGenCa)
  
# Superexpressos (principais pacientes)
summary(LungExpGenCaMax$SURVIVAL)
summary(SkinExpGenCaMax$SURVIVAL)

# Pouco expressos (principais pacientes)
summary(LungExpGenCaMin$SURVIVAL)
summary(SkinExpGenCaMin$SURVIVAL)

# Tratando a coluna SURVIVAL
LungExpGenCaMax <- mutate(LungExpGenCaMax, SURVIVAL=as.integer(SURVIVAL))
LungExpGenCaMin <- mutate(LungExpGenCaMin, SURVIVAL=as.integer(SURVIVAL))
SkinExpGenCaMax <- mutate(SkinExpGenCaMax, SURVIVAL=as.integer(SURVIVAL))
SkinExpGenCaMin <- mutate(SkinExpGenCaMin, SURVIVAL=as.integer(SURVIVAL))

LungExpGenCaMax <- arrange(LungExpGenCaMax, desc(SURVIVAL))
LungExpGenCaMin <- arrange(LungExpGenCaMin, desc(SURVIVAL))
SkinExpGenCaMax <- arrange(SkinExpGenCaMax, desc(SURVIVAL))
SkinExpGenCaMin <- arrange(SkinExpGenCaMin, desc(SURVIVAL))


## Separando os indivíduos com maior e menor expressão do gene CD274
# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, LungExpGenCa$Exp_Z_CD274 > median(LungExpGenCa$Exp_Z_CD274))
SkinExpGenCaMax <- filter(SkinExpGenCa, SkinExpGenCa$Exp_Z_CD274 > median(SkinExpGenCa$Exp_Z_CD274))
# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, LungExpGenCa$Exp_Z_CD274 < median(LungExpGenCa$Exp_Z_CD274))
SkinExpGenCaMin <- filter(SkinExpGenCa, SkinExpGenCa$Exp_Z_CD274 < median(SkinExpGenCa$Exp_Z_CD274))


#############################################################################
### Análise no LUNG

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)

### Análise no SKIN

SurvSkinMax <- select(SkinExpGenCaMax, SURVIVAL)
SurvSkinMax <- mutate(SurvSkinMax, Calda = "Superexpresso"); SurvSkinMax
SurvSkinMin <- select(SkinExpGenCaMin, SURVIVAL)
SurvSkinMin <- mutate(SurvSkinMin, Calda = "Subexpresso"); SurvSkinMin

## Survival Overall e Kaplan-Meier
SurvSkinMerge <- merge(SurvSkinMax, SurvSkinMin, all = TRUE); SurvSkinMerge
SurvSkinMerge$SURVIVAL <- as.double(SurvSkinMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvSkinMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvSkinMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4)
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvSkinMerge$Calda)))
plot(fit, col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)

####################################################################################################
####################################################################################################
#### LUNG
#### Selecionando TFs com p-valor conveniente para a análise.
#### Para isso, teremos que fazer o mesmo processo acima com todos os TFs.



### Quais genes?
genesL <- names(LungExpGenCa[6:dim(LungExpGenCa)[2]]); genesL

p.valueLung <- vector() 


####### Análise no LUNG --- GENE BACH1
i <- 1
genesL[i]
coluna <- LungExpGenCa$BACH1

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, LungExpGenCa$BACH1 < median(LungExpGenCa$BACH1))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE BCL3
i <- 2
genesL[i]
coluna <- LungExpGenCa$BCL3

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE CREM
i <- 3
genesL[i]
coluna <- LungExpGenCa$CREM

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE ETS1
i <- 4
genesL[i]
coluna <- LungExpGenCa$ETS1

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE ETS2
i <- 5
genesL[i]
coluna <- LungExpGenCa$ETS2

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE GATAD2A
i <- 6
genesL[i]
coluna <- LungExpGenCa$GATAD2A

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE HIVEP2
i <- 7
genesL[i]
coluna <- LungExpGenCa$HIVEP2

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE HMGA1
i <- 8
genesL[i]
coluna <- LungExpGenCa$HMGA1

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE HSF1
i <- 9
genesL[i]
coluna <- LungExpGenCa$HSF1

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE IRF1
i <- 10
genesL[i]
coluna <- LungExpGenCa$IRF1

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE NFKB2
i <- 11
genesL[i]
coluna <- LungExpGenCa$NFKB2

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE NFYA
i <- 12
genesL[i]
coluna <- LungExpGenCa$NFYA

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE PML
i <- 13
genesL[i]
coluna <- LungExpGenCa$PML

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE RARA
i <- 14
genesL[i]
coluna <- LungExpGenCa$RARA

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE REL
i <- 15
genesL[i]
coluna <- LungExpGenCa$REL

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE RELA
i <- 16
genesL[i]
coluna <- LungExpGenCa$RELA

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE RELB
i <- 17
genesL[i]
coluna <- LungExpGenCa$RELB

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE ETS1
i <- 19
genesL[i]
coluna <- LungExpGenCa$TEAD3

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE TEAD4
i <- 20
genesL[i]
coluna <- LungExpGenCa$TEAD4

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE PATZ1
i <- 21
genesL[i]
coluna <- LungExpGenCa$PATZ1

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE THAP11
i <- 22
genesL[i]
coluna <- LungExpGenCa$THAP11

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}


####### Análise no LUNG --- GENE ZNF554
i <- 23
genesL[i]
coluna <- LungExpGenCa$ZNF554

## Separando os indivíduos com maior e menor expressão

# Maior expressão
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)

# Menor expressão
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))

SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin

SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)

## Survival Overall e Kaplan-Meier
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
}







####################################################################################################
####################################################################################################
#### SKIN
#### Selecionando TFs com p-valor conveniente para a análise.
#### Para isso, teremos que fazer o mesmo processo acima com todos os TFs.
 
p.valueSkin <- vector()

### Quais genes?

genesS <- names(SkinExpGenCa[6:dim(SkinExpGenCa)[2]]); genesS

####### Análise no SKIN --- GENE RARG
j <- 1
genesS[j]
coluna <- SkinExpGenCa$RARG

### Análise no SKIN
## Separando os indivíduos com maior e menor expressão

# Maior expressão
SkinExpGenCaMax <- filter(SkinExpGenCa, coluna > median(coluna))

# Menor expressão
SkinExpGenCaMin <- filter(SkinExpGenCa, coluna < median(coluna))

SurvSkinMax <- select(SkinExpGenCaMax, SURVIVAL)
SurvSkinMax <- mutate(SurvSkinMax, Calda = "Superexpresso"); SurvSkinMax
SurvSkinMin <- select(SkinExpGenCaMin, SURVIVAL)
SurvSkinMin <- mutate(SurvSkinMin, Calda = "Subexpresso"); SurvSkinMin

## Survival Overall e Kaplan-Meier
SurvSkinMerge <- merge(SurvSkinMax, SurvSkinMin, all = TRUE); SurvSkinMerge
SurvSkinMerge$SURVIVAL <- as.double(SurvSkinMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvSkinMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvSkinMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4)
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvSkinMerge$Calda)))
plot(fit, main = as.name(genesS[j]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueSkin <- c(p.valueSkin, as.name(genesS[j])) 
}


####### Análise no SKIN --- GENE RXRA
j <- 2
genesS[j]
coluna <- SkinExpGenCa$RXRA

### Análise no SKIN
## Separando os indivíduos com maior e menor expressão

# Maior expressão
SkinExpGenCaMax <- filter(SkinExpGenCa, coluna > median(coluna))

# Menor expressão
SkinExpGenCaMin <- filter(SkinExpGenCa, coluna < median(coluna))

SurvSkinMax <- select(SkinExpGenCaMax, SURVIVAL)
SurvSkinMax <- mutate(SurvSkinMax, Calda = "Superexpresso"); SurvSkinMax
SurvSkinMin <- select(SkinExpGenCaMin, SURVIVAL)
SurvSkinMin <- mutate(SurvSkinMin, Calda = "Subexpresso"); SurvSkinMin

## Survival Overall e Kaplan-Meier
SurvSkinMerge <- merge(SurvSkinMax, SurvSkinMin, all = TRUE); SurvSkinMerge
SurvSkinMerge$SURVIVAL <- as.double(SurvSkinMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvSkinMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvSkinMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4)
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvSkinMerge$Calda)))
plot(fit, main = as.name(genesS[j]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueSkin <- c(p.valueSkin, as.name(genesS[j])) 
}


####### Análise no SKIN --- GENE TEF
j <- 3
genesS[j]
coluna <- SkinExpGenCa$TEF

### Análise no SKIN
## Separando os indivíduos com maior e menor expressão

# Maior expressão
SkinExpGenCaMax <- filter(SkinExpGenCa, coluna > median(coluna))

# Menor expressão
SkinExpGenCaMin <- filter(SkinExpGenCa, coluna < median(coluna))

SurvSkinMax <- select(SkinExpGenCaMax, SURVIVAL)
SurvSkinMax <- mutate(SurvSkinMax, Calda = "Superexpresso"); SurvSkinMax
SurvSkinMin <- select(SkinExpGenCaMin, SURVIVAL)
SurvSkinMin <- mutate(SurvSkinMin, Calda = "Subexpresso"); SurvSkinMin

## Survival Overall e Kaplan-Meier
SurvSkinMerge <- merge(SurvSkinMax, SurvSkinMin, all = TRUE); SurvSkinMerge
SurvSkinMerge$SURVIVAL <- as.double(SurvSkinMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvSkinMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvSkinMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4)
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvSkinMerge$Calda)))
plot(fit, main = as.name(genesS[j]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueSkin <- c(p.valueSkin, as.name(genesS[j])) 
}


####### Análise no SKIN --- GENE BATF3
j <- 4
genesS[j]
coluna <- SkinExpGenCa$BATF3

### Análise no SKIN
## Separando os indivíduos com maior e menor expressão

# Maior expressão
SkinExpGenCaMax <- filter(SkinExpGenCa, coluna > median(coluna))

# Menor expressão
SkinExpGenCaMin <- filter(SkinExpGenCa, coluna < median(coluna))

SurvSkinMax <- select(SkinExpGenCaMax, SURVIVAL)
SurvSkinMax <- mutate(SurvSkinMax, Calda = "Superexpresso"); SurvSkinMax
SurvSkinMin <- select(SkinExpGenCaMin, SURVIVAL)
SurvSkinMin <- mutate(SurvSkinMin, Calda = "Subexpresso"); SurvSkinMin

## Survival Overall e Kaplan-Meier
SurvSkinMerge <- merge(SurvSkinMax, SurvSkinMin, all = TRUE); SurvSkinMerge
SurvSkinMerge$SURVIVAL <- as.double(SurvSkinMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvSkinMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvSkinMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4)
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvSkinMerge$Calda)))
plot(fit, main = as.name(genesS[j]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueSkin <- c(p.valueSkin, as.name(genesS[j])) 
}


####### Análise no SKIN --- GENE IKZF1
j <- 5
genesS[j]
coluna <- SkinExpGenCa$IKZF1

### Análise no SKIN
## Separando os indivíduos com maior e menor expressão

# Maior expressão
SkinExpGenCaMax <- filter(SkinExpGenCa, coluna > median(coluna))

# Menor expressão
SkinExpGenCaMin <- filter(SkinExpGenCa, coluna < median(coluna))

SurvSkinMax <- select(SkinExpGenCaMax, SURVIVAL)
SurvSkinMax <- mutate(SurvSkinMax, Calda = "Superexpresso"); SurvSkinMax
SurvSkinMin <- select(SkinExpGenCaMin, SURVIVAL)
SurvSkinMin <- mutate(SurvSkinMin, Calda = "Subexpresso"); SurvSkinMin

## Survival Overall e Kaplan-Meier
SurvSkinMerge <- merge(SurvSkinMax, SurvSkinMin, all = TRUE); SurvSkinMerge
SurvSkinMerge$SURVIVAL <- as.double(SurvSkinMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvSkinMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvSkinMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4)
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvSkinMerge$Calda)))
plot(fit, main = as.name(genesS[j]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueSkin <- c(p.valueSkin, as.name(genesS[j])) 
}


####### Análise no SKIN --- GENE IRF7
j <- 6
genesS[j]
coluna <- SkinExpGenCa$IRF7

### Análise no SKIN
## Separando os indivíduos com maior e menor expressão

# Maior expressão
SkinExpGenCaMax <- filter(SkinExpGenCa, coluna > median(coluna))

# Menor expressão
SkinExpGenCaMin <- filter(SkinExpGenCa, coluna < median(coluna))

SurvSkinMax <- select(SkinExpGenCaMax, SURVIVAL)
SurvSkinMax <- mutate(SurvSkinMax, Calda = "Superexpresso"); SurvSkinMax
SurvSkinMin <- select(SkinExpGenCaMin, SURVIVAL)
SurvSkinMin <- mutate(SurvSkinMin, Calda = "Subexpresso"); SurvSkinMin

## Survival Overall e Kaplan-Meier
SurvSkinMerge <- merge(SurvSkinMax, SurvSkinMin, all = TRUE); SurvSkinMerge
SurvSkinMerge$SURVIVAL <- as.double(SurvSkinMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvSkinMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvSkinMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4)
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvSkinMerge$Calda)))
plot(fit, main = as.name(genesS[j]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueSkin <- c(p.valueSkin, as.name(genesS[j])) 
}


####### Análise no SKIN --- GENE PML
j <- 7
genesS[j]
coluna <- SkinExpGenCa$PML

### Análise no SKIN
## Separando os indivíduos com maior e menor expressão

# Maior expressão
SkinExpGenCaMax <- filter(SkinExpGenCa, coluna > median(coluna))

# Menor expressão
SkinExpGenCaMin <- filter(SkinExpGenCa, coluna < median(coluna))

SurvSkinMax <- select(SkinExpGenCaMax, SURVIVAL)
SurvSkinMax <- mutate(SurvSkinMax, Calda = "Superexpresso"); SurvSkinMax
SurvSkinMin <- select(SkinExpGenCaMin, SURVIVAL)
SurvSkinMin <- mutate(SurvSkinMin, Calda = "Subexpresso"); SurvSkinMin

## Survival Overall e Kaplan-Meier
SurvSkinMerge <- merge(SurvSkinMax, SurvSkinMin, all = TRUE); SurvSkinMerge
SurvSkinMerge$SURVIVAL <- as.double(SurvSkinMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvSkinMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvSkinMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4)
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvSkinMerge$Calda)))
plot(fit, main = as.name(genesS[j]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueSkin <- c(p.valueSkin, as.name(genesS[j])) 
}


####### Análise no SKIN --- GENE STAT1
j <- 8
genesS[j]
coluna <- SkinExpGenCa$STAT1

### Análise no SKIN
## Separando os indivíduos com maior e menor expressão

# Maior expressão
SkinExpGenCaMax <- filter(SkinExpGenCa, coluna > median(coluna))

# Menor expressão
SkinExpGenCaMin <- filter(SkinExpGenCa, coluna < median(coluna))

SurvSkinMax <- select(SkinExpGenCaMax, SURVIVAL)
SurvSkinMax <- mutate(SurvSkinMax, Calda = "Superexpresso"); SurvSkinMax
SurvSkinMin <- select(SkinExpGenCaMin, SURVIVAL)
SurvSkinMin <- mutate(SurvSkinMin, Calda = "Subexpresso"); SurvSkinMin

## Survival Overall e Kaplan-Meier
SurvSkinMerge <- merge(SurvSkinMax, SurvSkinMin, all = TRUE); SurvSkinMerge
SurvSkinMerge$SURVIVAL <- as.double(SurvSkinMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvSkinMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvSkinMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4)
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvSkinMerge$Calda)))
plot(fit, main = as.name(genesS[j]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueSkin <- c(p.valueSkin, as.name(genesS[j])) 
}


####### Análise no SKIN --- GENE STAT2
j <- 9
genesS[j]
coluna <- SkinExpGenCa$STAT2

### Análise no SKIN
## Separando os indivíduos com maior e menor expressão

# Maior expressão
SkinExpGenCaMax <- filter(SkinExpGenCa, coluna > median(coluna))

# Menor expressão
SkinExpGenCaMin <- filter(SkinExpGenCa, coluna < median(coluna))

SurvSkinMax <- select(SkinExpGenCaMax, SURVIVAL)
SurvSkinMax <- mutate(SurvSkinMax, Calda = "Superexpresso"); SurvSkinMax
SurvSkinMin <- select(SkinExpGenCaMin, SURVIVAL)
SurvSkinMin <- mutate(SurvSkinMin, Calda = "Subexpresso"); SurvSkinMin

## Survival Overall e Kaplan-Meier
SurvSkinMerge <- merge(SurvSkinMax, SurvSkinMin, all = TRUE); SurvSkinMerge
SurvSkinMerge$SURVIVAL <- as.double(SurvSkinMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvSkinMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvSkinMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4)
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvSkinMerge$Calda)))
plot(fit, main = as.name(genesS[j]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(4, 0.2, p)
if (p.val <= 0.05){
  p.valueSkin <- c(p.valueSkin, as.name(genesS[j])) 
}
#####################################################################################
################### CONCLUSÃO ESTATÍSTICA ###########################################
### Avaliando os vetores:

p.valueLung
p.valueSkin

# De acordo com a análise, somente o fator de transcrição BCL3 tem sua expressão gênica
# estatisticamente relacionada com a sobrevida dos pacientes.
# Visualisando o Kaplan-Meier do BCL3 novamente, temos:


####### Análise no LUNG --- GENE BCL3
i <- 2
genesL[i]
coluna <- LungExpGenCa$BCL3
LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
dim(LungExpGenCaMax)
LungExpGenCaMin <- filter(LungExpGenCa, coluna < median(coluna))
SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
# Plot do Kaplan-Meier
#jpeg("Kaplan-Meier - Gene BCL3.jpeg", pointsize = 20, width = 800, height = 800)
plot(fit, main = as.name(genesL[i]), col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
legend(65, 1, p)
#dev.off()


############################ fim ###############################################
  
############### LOOP ##################3


### Loop fail
for (i in seq(1,24)){
  genesL[i]
  coluna <- LungExpGenCa[i+6]
  ## Separando os indivíduos com maior e menor expressão
  
  # Maior expressão
  LungExpGenCaMax <- filter(LungExpGenCa, coluna > median(coluna))
  dim(LungExpGenCaMax)
  
  # Menor expressão
  LungExpGenCaMin <- filter(LungExpGenCa, LungExpGenCa$BACH1 < median(LungExpGenCa$BACH1))
  
  SurvLungMax <- select(LungExpGenCaMax, SURVIVAL)
  SurvLungMax <- mutate(SurvLungMax, Calda = "Superexpresso"); SurvLungMax
  SurvLungMin <- select(LungExpGenCaMin, SURVIVAL)
  SurvLungMin <- mutate(SurvLungMin, Calda = "Subexpresso"); SurvLungMin
  
  SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
  
  ## Survival Overall e Kaplan-Meier
  SurvLungMerge <- merge(SurvLungMax, SurvLungMin, all = TRUE); SurvLungMerge
  SurvLungMerge$SURVIVAL <- as.integer(SurvLungMerge$SURVIVAL)
  fit <- survfit(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
  diff <- survdiff(Surv(SURVIVAL) ~ Calda, data = SurvLungMerge)
  p.val <- round(1 - pchisq(diff$chisq, length(diff$n) - 1), digits = 4) # P-Valor do Kaplan-Meier
  p <- str_c("Valor-p: ",p.val)  # P-Valor do Kaplan-Meier
  Exp <- as.character(rev(unique(SurvLungMerge$Calda)))
  plot(fit, col=c(1:2), xlab="time(months)", ylab="Survival Probablity",lwd=1, conf.int = F)
  legend(23, 1, Exp, col=c(1:2), lwd = 0.3)
  legend(4, 0.2, p)
  if (p.val <= 0.05){
    p.valueLung <- c(p.valueLung, as.name(genesL[i])) 
  }
}
###
  
  

###########################################################################################
# NO GGplot


m <- length(fit)
df <- data.frame(time = c(rep(0, m), fit$time),
                 surv = c(rep(1, m), fit$surv),
                 group = c(names(fit$strata),
                           rep(names(fit$strata), fit$strata)))
ggplot(data = df) +
  geom_step(aes(x = time, y = surv, colour = as.factor(group))) +
  ylim(0, 1) +
  labs(colour = "Curvas", y = "Proporção de Sobreviventes", x = "Tempo")

