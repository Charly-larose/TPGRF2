######## TP GRF2, équipe 9 ######
## Émie Leclerc, Charliane Larose, Lorélie Gélinas

#### Installation des packages ####

liste <-c("magrittr", "derivmkts", "httr","ggplot2")

installation <- liste %in% installed.packages()
if(length(liste[!installation]) > 0) {
  install.packages(liste[!installation], repos = "https://cran.rstudio.com/")
}
lapply(liste, require, character.only = TRUE)

#### Importation des données ####
## importation taux banque du canada dans 5 dernières années
banque_can <- read.csv("taux_can.csv")


## importation données mensuelles historiques
data_histo <- read.csv("DonnéesTPGRF2(version2).csv")
data_histo <- data_histo[-1, 13]

## création de la fonction pour faire les arbres binomiaux ####
### pour voir le code source écrire nom de la foncton sans parenthèse
#View(binomplot)
arbre <- function (s, k, v, r, tt, d, nstep, putopt = FALSE, american = TRUE,
                   plotvalues = FALSE, plotarrows = FALSE, drawstrike = TRUE,
                   pointsize = 4, ylimval = c(0, 0), saveplot = FALSE, saveplotfn = "binomialplot.pdf",
                   crr = FALSE, jarrowrudd = FALSE, titles = TRUE, specifyupdn = FALSE,
                   up = 1.5, dn = 0.5, returnprice = FALSE, logy = FALSE)
{
    setylim <- ifelse((sum(ylimval^2) == 0), FALSE, TRUE)
    y <- binomopt(s, k, v, r, tt, d, nstep, american, putopt,
                  specifyupdn, crr, jarrowrudd, up, dn, returnparams = TRUE,
                  returntrees = TRUE)
    h <- tt/nstep
    for (i in c("up", "dn", "p")) assign(i, y$params[i])
    for (i in c("stree", "exertree", "oppricetree", "probtree")) assign(i,
                                                                        y[["i"]])
    nn <- 0:nstep
    payoffmult <- ifelse((putopt), -1, 1)
    stree <- y$stree
    exertree <- y$exertree
    oppricetree <- y$oppricetree
    probtree <- y$probtree
    bondtree <- y$bondtree
    deltatree <- y$deltatree
    plotcolor <- ifelse(exertree, "green3", "red")
    if (saveplot)
        pdf(saveplotfn)
    ylim_default <- c(0, max(stree) * 1.03)
    savepar <- par(no.readonly = TRUE)
    if (logy)
        par(ylog = TRUE)
    plot(rep(nn, nn + 1) * h, stree[stree > 0], ylim = ifelse(c(setylim,
                                                                setylim), ylimval, c(min(stree[stree > 1] - 0.95), max(stree) *
                                                                                         1.03)), col = plotcolor[stree > 0], pch = 21,
         cex = ifelse(stree[stree > 0], sqrt(probtree[stree > 0]) *
                          pointsize, 1), bg = plotcolor[stree > 0],
         xlab = ifelse(titles, "Période binomiale (année)", ""),
         ylab = ifelse(titles, "Prix de l'action ($)", ""),
         main = if (titles) paste(ifelse(putopt, "Option de vente", "Option d'achat"),
                                  ifelse(american, "américaine", "européenne")),
         log = ifelse(logy, "y", ""))

    if (titles)
        mtext(paste0("Stock = ", format(s, digits = 3), ",Prix d'exercice = ",
                     format(k, digits = 3), ", Temps = ", format(tt, digits = 4),
                     ifelse(tt == 1, " an,", " ans,"), " Prix = ",
                     format(oppricetree[1, 1], digits = 5)))
    if (drawstrike)
        abline(h = k, lty=2)
    yoffset <- ifelse(setylim, 0.03 * (ylimval[2] - ylimval[1]),
                      0.03 * max(stree))
    if (plotarrows) {
        for (i in 1:nstep) {
            for (j in 1:i) {
                arrows((i - 1) * h, stree[j, i], c(i, i) * h,
                       c(stree[j, i + 1], stree[j + 1, i + 1]), length = 0.06)
            }
        }
    }
    if (plotvalues) {
        for (i in 1:(nstep + 1)) {
            text((i - 1) * h, stree[1:i, i] + yoffset + 3, format(stree[1:i,
                                                                        i],
                                                                  digits = 3), cex = 0.8)
        }
        for (i in 1:(nstep)) {
            text((i - 1) * h, stree[1:i, i] + yoffset, format(deltatree[1:i,
                                                                        i],
                                                              digits = 3), cex = 0.8, col="deeppink")
        }

        for (i in 1:(nstep)) {
            text((i - 1) * h, stree[1:i, i] + yoffset-8, format(bondtree[1:i,
                                                                         i],
                                                                digits = 3), cex = 0.8, col="blue")
        }

        for (i in 1:(nstep)) {
            text((i - 1) * h, stree[1:i, i] + yoffset-11, format(oppricetree[1:i, i],
                                                                 digits = 3), cex = 0.8, col="purple")
        }
        legend("topleft", c("S", "Delta", "B", "Prix"), bty="n", lty=rep(1, 6),
               col=c("black", "deeppink", "blue", "purple"), cex=0.8)

    }
    if (saveplot)
        dev.off()
    if (returnprice)
        return(oppricetree[1, 1])
}


# test
# arbre(s = 100, k = 95, v = sig, r = r, tt = 1, d = 0, nstep = 4,
#       putopt = T, plotvalues = T, plotarrows = T, american = F,
#       returnprice = T)

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
### modèle à 4 périodes, illustrer les arbres

#### Dans l'arbre j'ai mis ligne pointillé pour prix d'exercice, est-ce qu'on enlève????
## Oui, je l'enleverais


## option de vente européenne, k = 95, prix initial de l'indice = 100r
r_log <- log(1+r) # transformation pour avoir un r composé continuellement

arbre(s = 100, k = 95, v = sig, r = r_log, tt = 1, d = 0, nstep = 4, american = FALSE,
      putopt = TRUE, plotvalues = TRUE, plotarrows = TRUE, returnprice = TRUE,
      drawstrike = TRUE)
## option d'achat américaine k = 110, prix initial de l'indice = 100
arbre(s = 100, k = 110, v = sig, r = r_log, tt = 1, d = 0, nstep = 4, american = TRUE,
      putopt = FALSE, plotvalues = TRUE, plotarrows = TRUE, returnprice = TRUE,
      drawstrike = TRUE)

### modèle à 52 périodes, seulement calculer les prix 
## option de vente européenne, k = 95, prix initial de l'indice = 100
binomopt(s = 100, k = 95, v = sig, r = r_log, tt = 1, d = 0, nstep = 52, american = F,
         putopt = T)
## option d'achat américaine k = 110, prix initial de l'indice = 100
binomopt(s = 100, k = 110, v = sig, r = r_log, tt = 1, d = 0, nstep = 52, american = T,
         putopt = F)




## option asiatiques
## option d'achat asiatique de type « option d'achat sur moyenne arithmétique »
# k = 110
# m aucune idée quoi mettre ????
arithavgpricecv(s = 100, k = 110, v = sig, r = r, t =1, d = 0, m = )
## option d'achat asiatique de type « option de vente à barrière désactivante »
# barrière = 105, k = 95
# Aucune idée comment faire, faudrait voir les notes


### Question 3 ###
## graphique pour illustrer la relation entre le prix d'exercice et le prix à
## payer pour les options d'achat et ventes européenne 

## CALL 
k <- seq(from = 50, to = 250, by = 5)
prix_achat <- sapply(k, function(k) binomopt(s = 100, k = k, v = sig, r = r_log, tt = 1, d = 0, nstep = 52, american = F,
                                        putopt = F, returntrees = T)$price)

data_graph <- data.frame(prix_exercice = k, prix_call = prix_achat)

graph_call <- ggplot(data_graph, aes(x = prix_exercice, y = prix_call)) + geom_line(linewidth = 1.2) +
  labs(
    title = "Prix d'un option d'achat européenne en fonction de différentes valeurs de prix d'exercice (L)",
    subtitle = "Le prix initial du sous-jacent est de 100 $ avec le modèle binomial à 52 périodes",
    x = "Prix d'exercice ($)",
    y = "Prix de l'option ($)",)+ theme_classic() + theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12)
    )

## PUT
k <- seq(from = 0, to = 200, by = 5)

prix_put <- sapply(k, function(k) binomopt(s = 100, k = k, v = sig, r = r_log, tt = 1, d = 0, nstep = 52, american = F,
                                             putopt = T, returntrees = T)$price)

data_graph <- data.frame(prix_exercice = k, prix_put = prix_put)

graph_put <- ggplot(data_graph, aes(x = prix_exercice, y = prix_put)) + geom_line(linewidth = 1.2) +
  labs(
    title = "Prix d'un option de vente européenne en fonction de différentes valeurs de prix d'exercice (k)",
    subtitle = "Le prix initial du sous-jacent est de 100 $ avec le modèle binomial à 52 périodes",
    x = "Prix d'exercice ($)",
    y = "Prix de l'option ($)",)+ theme_classic() + theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12)
    )






