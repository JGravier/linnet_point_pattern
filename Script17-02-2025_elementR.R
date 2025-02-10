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
# importer les données du projet SuDuCo
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
plot(x = rpoislpp(lambda = bijoutiers_infos$intensity, L = paris, nsim = 1),
     pch = 15, main = NULL)

# Écart entre observé et simulé (courbes)
# génération de 10 simulations
bijoutiers_10sim <- runiflpp(ex = bijoutiers_lpp, nsim = 10)

# calcul des distances entre les points simulés aléatoirement sur le réseau
list_simulated_dist <- list()

for (i in 1:length(bijoutiers_10sim)) {
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

library(ggplot2)
bijoutiers_ecart %>%
  ggplot(aes(x = dist_pi_p, group = n_sim)) +
  geom_density() +
  theme_bw()

bijoutiers_ecart %>%
  ggplot(aes(x = dist_pi_p, color = type)) +
  geom_density(linewidth = 1) +
  theme_bw()

# Écart à une répartition homogène : distribution de la longueur des plus courts chemins
# simuler 20 semis où distribution fonction longueur des voies
env <- envelope(LJ, linearK, correction="none", nsim=20)
plot(env)

# Écart à une répartition homogène
# au niveau de la ville (axe ouest-est)
# test de Berman
btB <- berman.test(bijoutiers_lpp, "x")
plot(btB)

# test de Kolmogorov-Smirnov
cdfB <- cdf.test(bijoutiers_lpp, "x")
plot(cdfB)

# distance au centre (test JG)

# Écart à une répartition homogène au niveau du tronçon
# épiciers davantage présents à proximité des intersections ?
alongE <- linfun(function(x,y,seg,tp) { tp }, domain(epiciers_lpp))
rhoalongE <- rhohat(epiciers_lpp, alongE)
plot(rhoalongE)

# hotspots
# calcul de la densité de points par unité de réseau
# l'option finespacing n'est pas nécessaire pour un réseau de petite taille
densite_bij <- density.lpp(x = bijoutiers_lpp, finespacing = FALSE, distance = "path")
densite_epic <- density.lpp(epiciers_lpp, finespacing = FALSE, distance = "path")

# option couleur (par défaut)
plot(densite_bij)

# option épaisseur des liens du réseau
plot(x = densite_epic, 
     main = "Densité de bijoutiers",
     style = "width", 
     adjust = 0.5)  # contrôler épaisseur max





