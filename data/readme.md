# Données contenues dans \data

## Réseau de rues de Jacoubet
### Référence
Les données spatiales "1836_jacoubet.shp" (format .shp et associés) sont issues de GeoHistoricalData. 2019. 
« “Verniquet map” and “Jacoubet Atlas” Paris street networks ». Harvard Dataverse. [https://doi.org/10.7910/DVN/CCESX4](https://doi.org/10.7910/DVN/CCESX4).

### Description
Tronçons des rues de Paris d'après le plan de Jacoubet de 1836.*

Chaque tronçon correspond à une polyligne composée de deux ou plus de points. 
Chaque intersection de rue implique la création d'un tronçon, ainsi une rue est composée de 1-N tronçons.

### Métadonnées

* id : identifiant du tronçon
* nom_entier : nom de la rue composée par le(s) tronçon(s)


## Bijoutiers et épiciers 
### Référence
Les données spatiales "jewellers_1839.gpkg" et "grocers_1839.gpkg" sont dérivées de GeoHistoricalData. 2023. 
« Annuaires historiques parisiens, 1798-1914. Extraction structurée et géolocalisée à l’adresse des 
listes nominatives par ordre alphabétique et par activité dans les volumes numérisés ». Geopackage. 
Dataset. NAKALA. [https://doi.org/10.34847/nkl.98eem49t](https://doi.org/10.34847/nkl.98eem49t).

### Description
Les bijoutiers (jewellers) et les épiciers (grocers) à Paris en 1839, d'après l'_Annuaire Général du Commerce_ de Paris de Lamy. 

Les points sont initialement géocodés automatiquement depuis les adresses des plans anciens de Paris (notamment le plan Jacoubet).

### Métadonnées

* fid : identifiant du point
* snap_dist : distance en mètres de la polyligne du réseau Jacoubet la plus proche sur laquelle est accroché le point
  (voir [maptools::snapPointsToLines()](https://www.rdocumentation.org/packages/maptools/versions/1.1-8/topics/snapPointsToLines))
* source : _Annuaire_ initial à partir duquel les données ont été créées, en l'occurrence l'_Annuaire Général du Commerce_ de Paris (Collection Didot),
  édité par Charles Lamy en 1839. Disponible sur Gallica : [https://gallica.bnf.fr/ark:/12148/bpt6k63243601/](https://gallica.bnf.fr/ark:/12148/bpt6k63243601/)
* source_annee : année de l'édition de l'_Annuaire_


## Pour en savoir plus sur la construction des données

Site du programme SoDUCo : [https://soduco.geohistoricaldata.org/](https://soduco.geohistoricaldata.org/)