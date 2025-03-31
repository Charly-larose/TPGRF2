#### Installation des packages ####



#### Importation des données ####
banque_can <- read.csv("taux_can.csv")

df <- read.csv("DonnéesTPGRF2(version2).csv")
df <- df[-1, c(1, 13)]
