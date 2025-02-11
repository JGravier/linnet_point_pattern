#L. Beauguitte et J. Gravier, 17 février 2025
#Analyse de semis de points sur un réseau spatial avec R
#Version provisoire à ne pas diffuser, merci !

#### Chargements des packages ####
library(sf) # package général de manipulation d'objets spatiaux
library(spatstat) # package général d'analyse spatiale de semis de points 2D
library(spatstat.linnet) # package particulier de la famille spatstat permettant d'étudier
# des semis de points sur un réseau planaire spatial
library(tibble) # package de gestion de tableaux (data frames)
library(tidyr) # package de manipulation de tableaux (forme et hiérarchie)
library(dplyr) # package de grammaire (verbe) de manipulation de tableaux
library(stringr) # package de manipulation de chaînes de caractères
library(ggplot2) # package de visualisation de données

#### Créer son propre jeu de données ####
##### Structure des objets du package spatstat.linnet #####
# Import des fichiers
# sommets et coordonnées (x,y)
mini_node <- read.table(file = "data/data_mini/mini_node.txt", 
                        header = TRUE, row.names = 1, sep = ",")

# extrémités des liens
mini_edge <- read.table("data/data_mini/mini_edge.txt", header = FALSE, sep = ",")

# transformation en objet ppp et matrice
mini_node_ppp <- ppp(x = mini_node$x, y = mini_node$y, c(0,10), c(0,10))
mini_edge_matrix <- as.matrix(mini_edge, ncol = 2)

# création d'un objet linnet : réseau planaire linéaire
mini <- linnet(vertices = mini_node_ppp, # sommets
               edges = mini_edge_matrix) # liens
class(mini)
plot(mini)

# premier semis de points sur l'affichage politique
semis_collage <- read.table("data/data_mini/mini_points1.txt", header = TRUE, sep = ",") 
semis_collage 
## coordonnées (x,y)
## nb : nombre d’affiches,
## pol : tendance politique (eg : extrême-gauche, fe : féministe, ed : extrême-droite)
## sti : présence d’autocollant (0 : non, 1 : oui)

# contrôle du typage des variables
str(semis_collage)

# typage de la variable sti en TRUE/FALSE
semis_collage$sti <- as.logical(semis_collage$sti)

# second semis de points sur les équipements publics
semis_equipement <- read.table("data/data_mini/mini_points2.txt", header = TRUE, sep = ",") 
## coordonnées (x,y)
## typologie (bus : arrêt de bus ; sub : station de métro ; ps : commissariat)


# transformation en semis de points sur réseau linéaire (objet lpp)
collage_lpp <- lpp(X = semis_collage, # semis de points
                   L = mini) # réseau planaire linéaire
equip_lpp <- lpp(X = semis_equipement, L = mini)


##### Visualisation d'un semis de points sur un réseau linéaire #####
# visualisation par défaut
plot(collage_lpp) # une visualisation par type de "marks" (attributs des points)

# choix de l'attribut à visualiser
plot(equip_lpp, which.marks = "pol")


#### Jeux de données spatiales (usuels) ####
##### Importer les données du projet SoDuCo #####
# import du réseau viaire parisien en 1836
paris <- st_read("data/1836_jacoubet.shp")

# import des semis de points en 1839
epiciers <- st_read(dsn = "data/grocers_1839.gpkg")
bijoutiers <- st_read(dsn = "data/jewellers_1839.gpkg")

# contrôle visuel des objets importés
plot(paris$geometry)
plot(epiciers$geom, pch = 15, col = "blue", add = TRUE)
plot(bijoutiers$geom, pch = 15, col = "red", bg = "red", add = TRUE)

##### Construction du réseau linéaire (linnet) et des semis de points (ppp) #####
# construction du réseau linéaire
paris <- as.psp(st_geometry(paris))
paris <- as.linnet(paris)

summary(paris)

# transformation en objet ppp
epiciers_ppp <- as.ppp(st_geometry(epiciers))
bijoutiers_ppp <- as.ppp(st_geometry(bijoutiers))


##### Propriétés de base des objets linnet et ppp #####
# nombre de sommets
nvertices(paris)

# longueur totale du réseau viaire (en mètres)
volume(paris)

# degré moyen
mean(vertexdegree(paris))

# Semis de points des épiciers
# nombre de points
npoints(epiciers_ppp)

# nombre de points par mètre
intensity(epiciers_ppp)
## appelée "intensité" dans le package
## revient à une densité par unité de base

# Semis de points des bijouteries
npoints(bijoutiers_ppp)
intensity(bijoutiers_ppp)

##### Calculs des distances entre les points #####
# distance géographique (euclidienne) entre paires de points
pairdist(epiciers_ppp)[1:5, 1:5] # matrice symétrique

# distance géographique (euclidienne) au plus proche voisin
nndist(epiciers_ppp)[1:5] # vecteur numérique, par défaut k = 1 (plus proche voisin)
## on peut spécifier le paramètre k pour obtenir les distances des kème voisins

# identifiant du plus proche voisin
nnwhich(epiciers_ppp)[1:5] # vecteur numérique

# intégrer le semis du point au réseau planaire (objet lpp)
epiciers_lpp <- lpp(X = epiciers_ppp, # semis
                    L = paris) # réseau
summary(epiciers_lpp)

# plus court chemin (géographique et non topologique) entre paires de points
pairdist.lpp(epiciers_lpp)[1:5, 1:5] # matrice symétrique
## NB: identique à l'utilisation de la fonction pairdist() sur un objet lpp

bijoutiers_lpp <- lpp(X = bijoutiers_ppp, L = paris) # création lpp
summary(bijoutiers_lpp)


# Écart à une simulation aléatoire (carte)
# visualisation d'une simulation
bijoutiers_infos <- summary(object = bijoutiers_lpp)
plot(x = rpoislpp(lambda = bijoutiers_infos$intensity, L = paris, nsim = 1),
     pch = 15, main = NULL)

# Écart entre observé et simulé (courbes)
# génération de 10 simulations
bijoutiers_10sim <- runiflpp(ex = bijoutiers_lpp, nsim = 10)

# calcul des distances entre les points simulés aléatoirement sur le réseau
list_simulated_dist <- list()

for (i in 1:length(bijoutiers_10sim)) { # boucle pour pédagogie
  # calculs des plus courts chemins sur réseau
  dist_pi_p <- spatstat.geom::pairdist(X = bijoutiers_10sim[[i]])
  dist_pi_p[upper.tri(x = dist_pi_p, diag = TRUE)] <- NA
  
  dist_pi_p <- dist_pi_p %>%
    tibble::as_tibble() %>%
    tibble::rowid_to_column(var = "Pi") %>%
    tidyr::pivot_longer(cols = -Pi, names_to = "P", values_to = "dist_pi_p") %>%
    dplyr::filter(!is.na(dist_pi_p)) %>%
    dplyr::mutate(P = stringr::str_replace_all(string = P, pattern = "V", replacement = "")) %>%
    dplyr::mutate(type = "simulation", n_sim = i)
  
  list_simulated_dist[[i]] <- dist_pi_p   
}

# création d'un unique tableau et non d'une liste de tableaux
table_simulated_dist <- do.call(what = "rbind", args = list_simulated_dist)
table_simulated_dist

# calcul des distances entre les points observés des bijoutiers en 1839
dist_bijoutiers <- spatstat.geom::pairdist(X = bijoutiers_lpp)
dist_bijoutiers[upper.tri(x = dist_bijoutiers, diag = TRUE)] <- NA

dist_bijoutiers <- dist_bijoutiers %>%
  as_tibble() %>%
  rowid_to_column(var = "Pi") %>%
  pivot_longer(cols = -Pi, names_to = "P", values_to = "dist_pi_p") %>%
  filter(!is.na(dist_pi_p)) %>%
  mutate(P = stringr::str_replace_all(string = P, pattern = "V", replacement = "")) %>%
  mutate(type = "observation", n_sim = NA)

bijoutiers_ecart <- table_simulated_dist %>%
  dplyr::bind_rows(dist_bijoutiers)

# Visualisation des écarts
bijoutiers_ecart %>%
  ggplot(aes(x = dist_pi_p, group = n_sim)) + # toutes les simulations individusalisées
  geom_density() +
  theme_bw()

bijoutiers_ecart %>%
  ggplot(aes(x = dist_pi_p, color = type)) +  # toutes les simulations groupées
  geom_density(linewidth = 1) +
  theme_bw()

# Écart à une répartition homogène : distribution de la longueur des plus courts chemins
# simuler 20 semis où distribution fonction longueur des voies
# CSR = complete spatial randomness
env <- envelope(bijoutiers_lpp, linearK, correction="none", nsim=20)
plot(env)

# Écart à une répartition homogène
# au niveau de la ville (axe ouest-est)
# test de Berman
btB <- berman.test(bijoutiers_lpp, "x")
plot(btB)

# test de Kolmogorov-Smirnov
cdfB <- cdf.test(bijoutiers_lpp, "x")
plot(cdfB)

# Écart à une répartition aléatoire en tenant compte de la distance au centre
# import des données
center <- st_read(dsn = "data/halles.gpkg") # Les Halles considéré comme centre économique
center_ppp <- as.ppp(st_geometry(center))
center_lpp <- lpp(X = center_ppp, L = paris)

# import fonction (adaptation de spatstat.linnet::distfun.lpp())
source(file = "local-functions.R")
f_dist2_center <- distfun.inverse.lpp(X = center_lpp) # construction de la fonction d'intensité

# répartition aléatoire, en fonction de l'inverse de la distance au centre 
# i.e. plus on s'éloigne, plus la probabilité qu'un point soit simulé sur un tronçon diminue
dist2_center_bijoutiers <- rlpp(n = nrow(bijoutiers), f = f_dist2_center, nsim = 1)
# dist2_center_bijoutiers <- rpoislpp(lambda = f_dist2_center, L = paris, nsim = 1)
plot(dist2_center_bijoutiers, pch = 15)

# Analyse des écarts
# génération de 10 simulations
bijoutiers_10sim <- rlpp(n = nrow(bijoutiers), f = f_dist2_center, nsim = 10)

# calcul des distances entre les points simulés entre les points simulés et le centre
list_simulated_dist <- list()

for (i in 1:length(bijoutiers_10sim)) {
  # calculs des plus courts chemins sur réseau
  dist_pi_p <- spatstat.geom::crossdist(X = bijoutiers_10sim[[i]], Y = center_lpp)
  dist_pi_p[upper.tri(x = dist_pi_p, diag = TRUE)] <- NA
  
  dist_pi_p <- dist_pi_p %>%
    tibble::as_tibble() %>%
    tibble::rowid_to_column(var = "Pi") %>%
    tidyr::pivot_longer(cols = -Pi, names_to = "P", values_to = "dist_pi_p") %>%
    dplyr::filter(!is.na(dist_pi_p)) %>%
    dplyr::mutate(P = stringr::str_replace_all(string = P, pattern = "V", replacement = "")) %>%
    dplyr::mutate(type = "simulation", n_sim = i)
  
  list_simulated_dist[[i]] <- dist_pi_p   
}

table_simulated_dist <- do.call(what = "rbind", args = list_simulated_dist)

# calcul des distances entre les points observés des bijoutiers et le centre
dist_bijoutiers <- spatstat.geom::crossdist(X = bijoutiers_lpp, Y = center_lpp)
dist_bijoutiers[upper.tri(x = dist_bijoutiers, diag = TRUE)] <- NA

dist_bijoutiers <- dist_bijoutiers %>%
  as_tibble() %>%
  rowid_to_column(var = "Pi") %>%
  pivot_longer(cols = -Pi, names_to = "P", values_to = "dist_pi_p") %>%
  filter(!is.na(dist_pi_p)) %>%
  mutate(P = stringr::str_replace_all(string = P, pattern = "V", replacement = "")) %>%
  mutate(type = "observation", n_sim = NA)

# Visualisation
bijoutiers_ecart <- table_simulated_dist %>%
  dplyr::bind_rows(dist_bijoutiers)

bijoutiers_ecart %>%
  ggplot(aes(x = dist_pi_p, color = type)) +
  geom_density(linewidth = 1) +
  theme_bw()

# Écart à une répartition homogène au niveau du tronçon
# épiciers davantage présents à proximité des intersections ?
alongE <- linfun(function(x,y,seg,tp) { tp }, domain(epiciers_lpp))
rhoalongE <- rhohat(epiciers_lpp, alongE)
plot(rhoalongE)

# hotspots
# calcul de la densité de points par unité de réseau
# l'option finespacing n'est pas nécessaire pour un réseau de petite taille
densite_bij <- density.lpp(x = bijoutiers_lpp, 
                           finespacing = FALSE, 
                           distance = "path")  # autre option : euclidian
densite_epic <- density.lpp(epiciers_lpp, 
                            finespacing = FALSE, 
                            distance = "path")

# option couleur (par défaut)
plot(densite_epic)

# option épaisseur des liens du réseau
plot(x = densite_bij, 
     main = "Densité de bijoutiers",
     style = "width", 
     adjust = 0.5)  # contrôler épaisseur max

# choix du lissage
densite_epi2 <- density.lpp(x = epiciers_lpp, 
                            sigma = 200,       # fonction de lissage en mètres
                            finespacing = FALSE, 
                            distance = "path") # distance plus court chemin
plot(densite_epi2)

# argument statistique pour le choix du lissage
# oblige à prendre en compte distance euclidienne

b <- bw.lppl(epiciers_lpp,
             distance = "e")
plot(b, main="Choix du seuil de lissage")

densite_epi3 <- density.lpp(x = epiciers_lpp, 
                            sigma = max(b), # fonction de lissage en mètres
                            finespacing = FALSE, 
                            distance = "path") # distance plus court chemin
plot(densite_epi3)





