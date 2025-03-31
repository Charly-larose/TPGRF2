######## TP GRF2, équipe 9 ######
## Émie Leclerc, Charliane Larose, Lorélie Gélinas

#### Installation des packages ####

liste <-c("magrittr", "derivmkts", "httr")

installation <- liste %in% installed.packages()
if(length(liste[!installation]) > 0) {
  install.packages(liste[!installation], repos = "https://cran.rstudio.com/")
}


#### Importation des données ####
## importation taux banque du canada dans 5 dernières années
banque_can <- read.csv("taux_can.csv")


## importation données mensuelles historiques
data_histo <- read.csv("DonnéesTPGRF2(version2).csv")
data_histo <- data_histo[-1, c(13)]

##### Question 1 ######
#### Approximation des paramètres arbre binomial ####

## moyenne taux banque du canada, mensuelle
r <- mean(banque_can[, 2])/100

h <- 1/12 # les données sont mensuelles

longeur_donnee <- length(data_histo)

donnee_ln <- numeric(longeur_donnee) # pour le remplir

donnee_ln[1] <- NA

for (i in 2:longeur_donnee){
  (donnee_ln[i] <- log(data_histo[i]/data_histo[i-1]))
}

sig <- sd(donnee_ln[-1])/sqrt(h)

u <- exp(r*h + sig*sqrt(h))
d <- exp(r*h - sig*sqrt(h))
p <- (exp(r*h)-d)/(u-d)


##### Question 2 ######



