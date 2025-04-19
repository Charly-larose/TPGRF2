##### ACT-2011 #####
##### Code Travail pratique #####
##### Équipe 9 ####
## Lorélie Gélinas, Charliane Larose, Émie Leclerc

#### Liste des paquetages et installation si nécessaire ####
liste <-c("magrittr", "derivmkts", "httr","ggplot2")

installation <- liste %in% installed.packages()
if(length(liste[!installation]) > 0) {
    install.packages(liste[!installation], repos = "https://cran.rstudio.com/")
}

lapply(liste, require, character.only = TRUE)

#### Fonction qui permet d'illustrer les arbres binomiaux ####
arbre <- function (s, k, v, r, tt, d, nstep, putopt = FALSE, american = TRUE,
                   plotvalues = FALSE, plotarrows = FALSE, drawstrike = TRUE,
                   pointsize = 4, ylimval = c(0, 0), saveplot = FALSE,
                   saveplotfn = "binomialplot.pdf",
                   crr = FALSE, jarrowrudd = FALSE, titles = TRUE,
                   specifyupdn = FALSE,
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
       setylim), ylimval,c(min(stree[stree > 1] - 0.95), max(stree) * 1.03)),
                                        col = plotcolor[stree > 0], pch = 21,
         cex = ifelse(stree[stree > 0], sqrt(probtree[stree > 0]) *
                          pointsize, 1), bg = plotcolor[stree > 0],
         xlab = ifelse(titles, "Période binomiale (année)", ""),
         ylab = ifelse(titles, "Prix de l'action ($)", ""),
         main = if (titles) paste(ifelse(putopt, "Option de vente",
                                         "Option d'achat"),
                                  ifelse(american, "américaine", "européenne")),
                                            log = ifelse(logy, "y", ""))

    if (titles)
        mtext(paste0("Prix initial = ", format(s, digits = 3),
                                ", Prix d'exercice = ",
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
            text((i - 1) * h, stree[1:i, i] + yoffset + 3, format(stree[1:i, i],
                                                      digits = 3), cex = 0.8)
        }
        for (i in 1:(nstep)) {
            text((i - 1) * h, stree[1:i, i] + yoffset, format(deltatree[1:i, i],
                                        digits = 3), cex = 0.8, col="deeppink")
        }

        for (i in 1:(nstep)) {
            text((i - 1) * h, stree[1:i, i] + yoffset-8, format(bondtree[1:i,
                                                                         i],
                                            digits = 3), cex = 0.8, col="blue")
        }

        for (i in 1:(nstep)) {
            text((i - 1) * h, stree[1:i, i] + yoffset-11, format(oppricetree[1:i,
                                                                             i],
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


#### 1. Approximation des paramètres ####
## importation des données pour approximer le taux sans risque
banque_can <- read.csv("taux_can.csv")

## importation des données mensuelles historiques
data_histo <- read.csv("DonnéesTPGRF2(version2).csv")
data_histo <- data_histo[-1, 13]

## moyenne taux banque du canada, base effectif annuel
r <- mean(banque_can[, 2])/100
## convertir le taux effectif annuel en taux composé continuellement
r_log <- log(1 + r)

h <- 1/12 # les données historiques sont mensuelles

## approximation de la voltilité
# création d'un vecteur
longeur_donnee <- length(data_histo)

donnee_ln <- numeric(longeur_donnee) # pour le remplir

donnee_ln[1] <- NA

# remplir le vecteur
for (i in 2:longeur_donnee){
    (donnee_ln[i] <- log(data_histo[i]/data_histo[i-1]))
}

# appliquation de la formule pour obtenir l'approximation de la volatilité
sig <- sd(donnee_ln[-1])/sqrt(h)

## Paramètre : u, d, et p pour l'arbre binomial à 4 périodes
u_quatreper <- exp((r_log*1/4)+sig*sqrt(1/4))
d_quatreper <- exp((r_log*1/4)-sig*sqrt(1/4))
p_quatreper <- (exp(r_log*1/4)-d_quatreper)/(u_quatreper-d_quatreper)

## Paramètre : u, d, et p pour l'arbre binomial à 52 périodes
u_52per <- exp((r_log*1/52)+sig*sqrt(1/52))
d_52per <- exp((r_log*1/52)-sig*sqrt(1/52))
p_52per <- (exp(r_log*1/52)-d_52per)/(u_52per-d_52per)


#### 2. Arbres binomiaux ####
### Arbres binomiaux à 4 périodes ###
## Option de vente européenne
# Illustration de l'arbre de l'option de vente européenne
arbre(s = 100, k = 95, v = sig, r = r_log, tt = 1, d = 0, nstep = 4,
      american = FALSE, putopt = TRUE, plotvalues = TRUE, plotarrows = TRUE,
                    returnprice = TRUE,drawstrike = FALSE)
# Mettre le prix dans un objet pour pouvoir l'utiliser dans le texte
put95_4p = binomopt(s = 100, k = 95, v = sig, r = r_log, tt = 1, d = 0,
                                nstep = 4, american = F, putopt = T)

## Option d'achat européenne
# Illustration de l'arbre de l'option d'achat européenne
arbre(s = 100, k = 110, v = sig, r = r_log, tt = 1, d = 0, nstep = 4,
      american = FALSE, putopt = FALSE, plotvalues = TRUE, plotarrows = TRUE,
                        returnprice = TRUE, drawstrike = FALSE)
# Mettre le prix dans un objet pour pouvoir l'utiliser dans le texte
call110_4p = binomopt(s = 100, k = 110, v = sig, r = r_log, tt = 1, d = 0,
                      nstep = 4, american = F, putopt = F)


### Arbres binomiaux à 52 périodes ###
## Prix option de vente européenne
put95_52p = binomopt(s = 100, k = 95, v = sig, r = r_log, tt = 1, d = 0,
                     nstep = 52, american = F, putopt = T)

## Prix option d'achat européenne
call110_52p = binomopt(s = 100, k = 110, v = sig, r = r_log, tt = 1, d = 0,
                       nstep = 52, american = F, putopt = F)


## Prix option d'achat asiatique
set.seed(20030811)
call_asia <- arithavgpricecv(100, 110, sig, r_log, 1, 0, m=52,
                                                    numsim=1000000)[[1]]

## Prix option de vente à barrière désactivante
put_asia <-putupout(100, 95, sig, r_log, 1, 0, 105)


#### 3. Relation du prix de l'option et du prix d'exercice ####
## Option d'achat
# prix d'exercice pour lequel on veut calculer le prix de l'option
k <- seq(from = 50, to = 250, by = 5)
# fonction pour calculer le prix de l'option
prix_achat <- sapply(k, function(k) binomopt(s = 100, k = k, v = sig, r = r_log,
                                        tt = 1, d = 0, nstep = 52, american = F,
                                             putopt = F, returntrees = T)$price)
# code pour faire le graphique
data_graph <- data.frame(prix_exercice = k, prix_call = prix_achat)

graph_call <- ggplot(data_graph, aes(x = prix_exercice, y = prix_call)) +
    geom_line() +
    labs(
        x = "Prix d'exercice ($)",
        y = "Prix de l'option ($)"
    ) +
    theme(
        plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12)
    ) + theme_classic()
# afficher le graphique
graph_call

## Option de vente
# prix d'exercice pour lequel on veut calculer le prix de l'option
k <- seq(from = 0, to = 200, by = 5)
# fonction pour calculer le prix de l'option
prix_put <- sapply(k, function(k) binomopt(s = 100, k = k, v = sig, r = r_log,
                                        tt = 1, d = 0, nstep = 52, american = F,
                                           putopt = T, returntrees = T)$price)
# code pour faire le graphique
data_graph <- data.frame(prix_exercice = k, prix_put = prix_put)

graph_put <- ggplot(data_graph, aes(x = prix_exercice, y = prix_put)) +
    geom_line() +
    labs(
        x = "Prix d'exercice ($)",
        y = "Prix de l'option ($)"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12))
# afficher le graphique
graph_put


#### 4. Abre binomial de l'option de vente américaine ####
## Illustration de l'arbre
arbre(s = 100, k = 95, v = sig, r = r_log, tt = 1, d = 0, nstep = 4,
      american = TRUE, putopt = TRUE, plotvalues = TRUE, plotarrows = TRUE,
                    returnprice = TRUE, drawstrike = FALSE)
## prix de l'option afin de pouvoir l'utiliser dans le rapport
put95_4p_ame = binomopt(s = 100, k = 95, v = sig, r = r_log, tt = 1, d = 0,
                        nstep = 4, american = T, putopt = T)

#### 5. Modèle de Black-Scholes ####
## Option d'achat
call_BS <- function(s, k, sig, r_log, t, div)
{
    d1 <- (log(s/k) + (r_log - div + sig^2/2)*t)/(sig*sqrt(t))
    d2 <- d1 - sig*sqrt(t)
    s*exp(-div*t)*pnorm(d1) - k*exp(-r_log*t)*pnorm(d2)
}
# Prix
call_BS110 <- call_BS(100, 110, sig, r_log, 1, 0)

## Option de vente
put_BS <- function(s, k, sig, r_log, t, div)
{
    d1 <- (log(s/k) + (r_log - div + sig^2/2)*t)/(sig*sqrt(t))
    d2 <- d1 - sig*sqrt(t)
    k*exp(-r_log*t)*pnorm(-d2) - s*exp(-div*t)*pnorm(-d1)
}
# Prix
put_BS95 <- put_BS(100, 95, sig, r_log, 1, 0)


## Comparaison entre le modèle de l'arbre binomial et le modèle de Black-Scholes
# code pour montrer la converge du modèle de l'arbre binomial et le modèle de
# Black-Scholes pour l'option d'achat
s     <- 100
k     <- 110
tt    <- 1
d     <- 0


bs_call <- call_BS(s = s, k = k, sig = sig, r_log = r_log, t = tt, div = d)

# Boucle sur n pas pour le binomial
n_max       <- 100
binom_prices_call <- sapply(1:n_max, function(n) {
    arbre(
        s           = s,
        k           = k,
        v           = sig,
        r           = r_log,
        tt          = tt,
        d           = d,
        nstep       = n,
        putopt      = FALSE,
        american    = FALSE,
        returnprice = TRUE
    )
})


df <- data.frame(
    nstep    = 1:n_max,
    Binomial = binom_prices_call
)


convergence_call <- ggplot(df, aes(x = nstep, y = Binomial)) +
    geom_line(color = "forestgreen", size = 1) +
    geom_hline(yintercept = bs_call,
               linetype    = "dashed",
               color       = "red",
               size        = 1) +
    annotate("text",
             x     = round(n_max * 0.8),
             y     = bs_call + 0.05,
             label = "Black–Scholes",
             color = "red") +
    labs(
        x     = "Nombre de périodes",
        y     = "Prix de l'option"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
convergence_call


# code pour montrer la converge du modèle de l'arbre binomial et le modèle de
# Black-Scholes pour l'option d'achat
s     <- 100
k     <- 95
tt    <- 1
d     <- 0


bs_put <- put_BS(s = s, k = k, sig = sig, r_log = r_log, t = tt, div = d)

# Boucle sur n pas pour le binomial
n_max       <- 100
binom_prices <- sapply(1:n_max, function(n) {
    arbre(
        s           = s,
        k           = k,
        v           = sig,
        r           = r_log,
        tt          = tt,
        d           = d,
        nstep       = n,
        putopt      = TRUE,
        american    = FALSE,
        returnprice = TRUE
    )
})


df <- data.frame(
    nstep    = 1:n_max,
    Binomial = binom_prices
)


convergence <- ggplot(df, aes(x = nstep, y = Binomial)) +
    geom_line(color = "forestgreen", size = 1) +
    geom_hline(yintercept = bs_put,
               linetype    = "dashed",
               color       = "red",
               size        = 1) +
    annotate("text",
             x     = round(n_max * 0.8),
             y     = bs_put + 0.05,
             label = "Black–Scholes",
             color = "red") +
    labs(x     = "Nombre de périodes",
         y     = "Prix de l'option"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
convergence


# coût de l'option d'achat européenne en fonction du prix d'exercice
# avec le modèle de Black-Scholes
k <- seq(from = 50, to = 250, by = 5)
prix_achat_BS <- sapply(k, function(k) call_BS(s = s, k = k, sig = sig,
                                               r_log = r_log, t = tt, div = d))

data_graph_BS <- data.frame(prix_exercice = k, prix_call_BS = prix_achat_BS)

graph_call_BS <- ggplot(data_graph_BS, aes(x = prix_exercice, y = prix_call_BS))
    + geom_line() +
    labs(
        x = "Prix d'exercice ($)",
        y = "Prix de l'option ($)"
    ) +
    theme(
        plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12)
    ) + theme_classic()
graph_call_BS

# coût de l'option de vente européenne en fonction du prix d'exercice
# avec le modèle de Black-Scholes
k <- seq(from = 0, to = 200, by = 5)
prix_put_BS <- sapply(k, function(k) put_BS(s = s, k = k, sig = sig,
                                            r_log = r_log, t = tt, div = d))

data_graph_BS <- data.frame(prix_exercice = k, prix_put_BS = prix_put_BS)

graph_put_BS <- ggplot(data_graph_BS, aes(x = prix_exercice, y = prix_put_BS)) +
    geom_line() +
    labs(
        x = "Prix d'exercice ($)",
        y = "Prix de l'option ($)"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12))
graph_put_BS

