#L. Beauguitte et J. Gravier, 17 février 2025
#Analyse de semis de points sur un réseau spatial avec R
#Version provisoire à ne pas diffuser, merci !

# chargements des packages
library(sf)
library(spatstat)
library(tibble)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

# créer son propre jeu de données
# import des fichiers
# sommets et coordonnées
mini_node <- read.table(file = "data/data_mini/mini_node.txt", header = TRUE, row.names = 1, sep = ",")

# extrémités des liens
mini_edge <- read.table("data/data_mini/mini_edge.txt", header = FALSE, sep = ",")

# transformation en objet ppp et matrice
mini_node_ppp <- ppp(x = mini_node$x, y = mini_node$y, c(0,10), c(0,10))
mini_edge_matrix <- as.matrix(mini_edge, ncol = 2)

mini <- linnet(vertices = mini_node_ppp, edges = mini_edge_matrix)
class(mini)
plot(mini)

# premier semis sur l'affichage politique
semis_collage <- read.table("data/data_mini/mini_points1.txt", header = TRUE, sep = ",") 

# contrôle du typage des variables
str(semis_collage)

# typage de la variable sti
semis_collage$sti <- as.logical(semis_collage$sti)

# second semis sur les équipements publics
semis_equipement <- read.table("data/data_mini/mini_points2.txt", header = TRUE, sep = ",") 

# transformation en semis de points sur réseau linéaire
collage_lpp <- lpp(X = semis_collage, L = mini)
equip_lpp <- lpp(X = semis_equipement, L = mini)

# visualisation par défaut
plot(collage_lpp)

# choix de l'attribut à visualiser
plot(equip_lpp, which.marks = "pol")

#######################################
# importer les données du projet SoDuCo
# import du réseau viaire parisien en 1836
paris <- st_read("data/1836_jacoubet.shp")

# import des semis de points
epiciers <- st_read(dsn = "data/grocers_1839.gpkg")
bijoutiers <- st_read(dsn = "data/jewellers_1839.gpkg")

# contrôle visuel des objets importés
plot(paris$geometry)
plot(epiciers$geom, pch = 15, col = "blue", add = TRUE)
plot(bijoutiers$geom, pch = 15, col = "red", bg = "red", add = TRUE)

# transformation en objet ppp
epiciers_ppp <- as.ppp(st_geometry(epiciers))
bijoutiers_ppp <- as.ppp(st_geometry(bijoutiers))

paris <- as.psp(st_geometry(paris))
paris <- as.linnet(paris)

summary(paris)

# Propriétés de base
# nombre de sommets
nvertices(paris)

# longueur totale du réseau viaire (en mètres)
volume(paris)

# degré moyen
mean(vertexdegree(paris))

# semis des épiciers
# nombre de points et nombre de points par mètre
npoints(epiciers_ppp)
intensity(epiciers_ppp)

# plus court chemin (géographique et non topologique) entre paires de points
pairdist(epiciers_ppp)[1:5, 1:5]

# distance géographique au plus proche voisin
nndist(epiciers_ppp)[1:5] 

# identifiant du plus proche voisin
nnwhich(epiciers_ppp)[1:5] 

# semis des bijouteries
npoints(bijoutiers_ppp)
intensity(bijoutiers_ppp)

# intégrer semis du point au réseau planaire
epiciers_lpp <- lpp(X = epiciers_ppp, L = paris)
summary(epiciers_lpp)

bijoutiers_lpp <- lpp(X = bijoutiers_ppp, L = paris)
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

# Écart à une répartition aléatoire tenant compte de la distance au centre
# import des données
center <- st_read(dsn = "data/halles.gpkg") # Les Halles considéré comme centre économique
center_ppp <- as.ppp(st_geometry(center))
center_lpp <- lpp(X = center_ppp, L = paris)

# import fonction (adaptation de spatstat.linnet::distfun.lpp())
source(file = "local-functions.R")
f_dist2_center <- distfun.inverse.lpp(X = center_lpp) # construction fonction

# répartition aléatoire, fonction de l'inverse de la distance au centre 
# (i.e. plus on s'éloigne, plus la probabilité qu'un point soit simulé sur un tronçon diminue)
dist2_center_bijoutiers <- rlpp(n = nrow(bijoutiers), f = f_dist2_center, nsim = 1)
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





