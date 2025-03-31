######## TP GRF2, équipe 9 ######
## Émie Leclerc, Charliane Larose, Lorélie Gélinas

#### Installation des packages ####
REQUIREDPACKAGES <- c("magrittr", "derivmkts", "httr")
## Install (if not installed) and load necessary packages.
for (package in REQUIREDPACKAGES) {
    if (!package %in% installed.packages()[, 1])
        install.packages(package)
    eval(bquote(library(.(package))))
}

library(derivmkts)

#### Importation des données ####
## importation taux banque du canada dans 5 dernières années
banque_can <- read.csv("taux_can.csv")
## moyenne taux banque du canada, mensuelle
moy_taux_can <- mean(banque_can[, 2])/100

## importation données mensuelles historiques
data_histo <- read.csv("DonnéesTPGRF2(version2).csv")
data_histo <- data_histo[-1, c(1, 13)]

#### Approximation des paramètres arbre binomial ####








