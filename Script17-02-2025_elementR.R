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
nndist(epiciers_ppp)[1:5] # vecteur numérique ; par défaut k = 1 (plus proche voisin)
## on peut spécifier le paramètre k pour obtenir les distances des kème voisins

# identifiant du plus proche voisin
nnwhich(epiciers_ppp)[1:5] # vecteur numérique

# intégrer le semis du point au réseau planaire (objet lpp)
epiciers_lpp <- lpp(X = epiciers_ppp, # semis
                    L = paris) # réseau
summary(epiciers_lpp)

# distance au plus court chemin (géographique et non topologique) entre paires de points
pairdist.lpp(epiciers_lpp)[1:5, 1:5] # matrice symétrique
## pairdist.lpp() identique à l'utilisation de pairdist() sur un objet lpp

# distance géographique (au plus court chemin) au plus proche voisin
nndist.lpp(epiciers_lpp, k = 2)[1:5] # vecteur numérique
## nndist.lpp() identique à l'utilisation de nndist() sur un objet lpp

# intégrer le semis du point au réseau planaire des bijoutiers
bijoutiers_lpp <- lpp(X = bijoutiers_ppp, L = paris) # création lpp


#### Écart à une répartition spatiale aléatoire ####
##### Visualisation cartographique d'une simulation aléatoire homogène #####
# Simulation homogène de Poisson d'un semis de points sur réseau selon l'intensité moyenne
bijoutiers_infos <- summary(object = bijoutiers_lpp) # pour récupérer l'intensité moyenne
plot(x = rpoislpp(lambda = bijoutiers_infos$intensity, # simulation de Poisson
                  # lambda est une constante car simulation homogène
                  L = paris, # réseau linéaire
                  nsim = 1), # nombre de simulation, par défaut nsim = 1
     pch = 15, # type de points cartographiés (carrés)
     main = NULL) # titre de la carte (ici, aucun)


##### Distribution des distances des plus courts chemins (observés vs simulés): courbes #####
# génération de 10 simulations
bijoutiers_10sim <- rpoislpp(ex = bijoutiers_lpp, # paramètre permettant de ne pas spécifier lambda et L
                             # alors, lambda = intensité moyenne du lpp
                             # L est le réseau de lpp
                             nsim = 10)
bijoutiers_10sim # une liste de 10 lpp simulés

# calcul des distances entre les points simulés aléatoirement sur le réseau
list_simulated_dist <- list()

for (i in 1:length(bijoutiers_10sim)) { # boucle pour pédagogie
  # calculs des distances des plus courts chemins sur réseau
  dist_pi_p <- pairdist.lpp(X = bijoutiers_10sim[[i]]) # matrice symétrique
  dist_pi_p[upper.tri(x = dist_pi_p, diag = TRUE)] <- NA # transformation de la partie
  # triangulaire haute de la matrice et de la diagonale en NA
  
  dist_pi_p <- dist_pi_p %>%
    tibble::as_tibble() %>% # matrice en tableau
    tibble::rowid_to_column(var = "Pi") %>% # ajout d'identifiant dans une colonne Pi
    tidyr::pivot_longer(cols = -Pi, # transformation du tableau en format long
                        # toutes les colonnes sauf Pi
                        names_to = "P", # variable contenant les noms des 1 à N colonnes initiales de la matrice
                        values_to = "dist_pi_p") %>% # variable contenant les valeurs de distances
    dplyr::filter(!is.na(dist_pi_p)) %>% # suppression de toutes les lignes contenant des NA
    dplyr::mutate(P = stringr::str_replace_all(string = P, # suppression de tous les "V" des noms des colonnes initiales de la matrice
                                               pattern = "V", 
                                               replacement = "")) %>%
    dplyr::mutate(type = "simulation", # nouvelle variable type
                  n_sim = i) # variable contenant le numéro de la simulation
  
  # liste de tableaux longs de distance au plus court chemin
  list_simulated_dist[[i]] <- dist_pi_p # 1 tableau pour 1 simulation
}

# création d'un unique tableau et non d'une liste de tableaux
table_simulated_dist <- do.call(what = "rbind", args = list_simulated_dist)
table_simulated_dist

# calcul des distances entre les points observés des bijoutiers en 1839
dist_bijoutiers <- pairdist.lpp(X = bijoutiers_lpp)
dist_bijoutiers[upper.tri(x = dist_bijoutiers, diag = TRUE)] <- NA

dist_bijoutiers <- dist_bijoutiers %>%
  as_tibble() %>%
  rowid_to_column(var = "Pi") %>%
  pivot_longer(cols = -Pi, names_to = "P", values_to = "dist_pi_p") %>%
  filter(!is.na(dist_pi_p)) %>%
  mutate(P = stringr::str_replace_all(string = P, pattern = "V", replacement = "")) %>%
  mutate(type = "observation", n_sim = NA)

# construction d'un tableau regroupant les simulations et observations
bijoutiers_ecart <- table_simulated_dist %>% # tableau des distances calculées pour les simulations
  dplyr::bind_rows(dist_bijoutiers) # ajout à la suite des lignes du tableau des observations
bijoutiers_ecart

# Visualisation graphique des écarts
bijoutiers_ecart %>%
  ggplot(aes(x = dist_pi_p, # x : représentation des distances
             group = n_sim)) + # toutes les simulations individusalisées
  geom_density() + # courbe de densité (KDE)
  theme_bw() # thème visuel fond clair

bijoutiers_ecart %>%
  ggplot(aes(x = dist_pi_p, 
             color = type)) +  # toutes les simulations sont regroupées en 1 courbe
  geom_density(linewidth = 1) + # plus forte épaisseur des traits
  xlab("distance en mètres") + # titre de l'axe des x
  ylab("KDE") + # titre de l'axe des y
  theme_bw()


##### Distribution des distances: Fonction K #####
# depuis chacun des points du semis au sein d'enveloppes de distance sur réseau (noté r)

# simuler 20 semis aléatoires inhomogènes
# où la distribution des points est fonction de la longueur des tronçons
# i.e.: la probabilité de titer un tronçon Ls sur L est fonction de sa longueur dans L
# le point est généré selon une probabilité uniforme le long de Ls
# c'est ce qu'on appelle une CSR = complete spatial randomness
env <- envelope.lpp(Y = bijoutiers_lpp, # objet lpp observé
                    fun = linearK, # fonction calculée pour chaque semis (simulé et observé)
                    correction = "none", # pas de correction de la géométrie du réseau
                    nsim = 20)
plot(env)

##### Intégrer des hypothèses spatiales complémentaires dans les analyses: échelle urbaine #####
# Écart à une répartition homogène (Poisson) à l'échelle urbaine selon un axe ouest-est
# test de Berman
bt_berman <- berman.test(X = bijoutiers_lpp, covariate = "x")
plot(bt_berman)

# test de Kolmogorov-Smirnov
ks <- cdf.test(bijoutiers_lpp, "x")
plot(ks)

# Écart à une répartition aléatoire inhomogène qui tienne compte de la distance au centre
# import des données
center <- st_read(dsn = "data/halles.gpkg") # Les Halles considérées comme centre économique
center_ppp <- as.ppp(st_geometry(center)) # semis de points
center_lpp <- lpp(X = center_ppp, L = paris) # semis de points sur le réseau

# calcul et visualisation de la distance au centre
plot(distfun.lpp(X = center_lpp))

# import fonction
source(file = "local-functions.R") # adaptation de spatstat.linnet::distfun.lpp()
f_dist2_center <- distfun.inverse.lpp(X = center_lpp) # création de la fonction d'intensité
f_dist2_center

# répartition aléatoire inhomogène selon l'inverse de la distance au centre 
# i.e. plus on s'éloigne, plus la probabilité qu'un point soit simulé sur un tronçon diminue
# simulation aléatoire de Poisson
dist2_center_bijoutiers <- rpoislpp(lambda = f_dist2_center, # fonction d'intensité
                                    L = paris, # réseau
                                    lmax = 1,
                                    nsim = 1)
# visualisation
plot(dist2_center_bijoutiers, pch = 15)

# autre simulation aléatoire
## NB: utilisé ici pour une question de temps de calcul
dist2_center_bijoutiers <- rlpp(n = nrow(bijoutiers), # nombre de points que l'on veut générer
                                f = f_dist2_center, # fonction d'intensité (densité de probabilité)
                                nsim = 1)
plot(dist2_center_bijoutiers, pch = 15)

# Analyse des écarts (observés vs simulés)
# génération de 10 simulations
bijoutiers_10sim <- rlpp(n = nrow(bijoutiers), f = f_dist2_center, nsim = 10)

# calcul des distances entre les points simulés et le centre
list_simulated_dist <- list()

for (i in 1:length(bijoutiers_10sim)) {
  dist_pi_p <- crossdist.lpp(X = bijoutiers_10sim[[i]], Y = center_lpp) # distances au plus court chemin entre deux lpp
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
dist_bijoutiers <- crossdist.lpp(X = bijoutiers_lpp, Y = center_lpp)
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
  xlab("distance en mètres") +
  ylab("KDE") +
  theme_bw()

##### Intégrer des hypothèses complémentaires: échelle des tronçons #####
# Écart à une répartition homogène au niveau du tronçon
# épiciers davantage présents à proximité des intersections ?
along_edges <- linfun(function(x,y,seg,tp) { tp }, domain(epiciers_lpp))
rhoalong_edges <- rhohat(object = epiciers_lpp, covariate = along_edges)
plot(rhoalong_edges)

#### Concentrations spatiales ####
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





