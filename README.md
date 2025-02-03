# Analyse de semis de points sur un réseau spatial avec R

_Document en préparation pour séance [ElementR](https://elementr.gitpages.huma-num.fr/website/apropos.html) 17 février 2025_

**Résumé** : L’analyse d’un semis de points sur un réseau désigne un ensemble de méthodes statistiques permettant de caractériser des événements ponctuels dans le temps et dans l’espace prenant place sur un réseau spatial planaire. Deux des sujets les plus traités dans la bibliographie sont les accidents de la circulation et les actes criminels commis dans l’espace public. Il est possible d’imaginer étudier d’autres thématiques plus ou moins ponctuelles dans le temps, qu’il s’agisse de l’offre commerciale dans un espace donné, du collage militant ou publicitaire, de la présence d’équipements dans l’espace public (bancs, toilettes), etc.

Trois grands types de question sont généralement posés à ces données :
   - Les points sont-ils significativement proches (ou éloignés) les uns des autres (étude du voisinage) ?
   - Existe-t-il des lieux où la concentration des points est notable (hot spots) ?
   - Quel modèle statistique est susceptible d’expliquer la géographie de ce semis de points ?

Ces trois questions seront abordées durant la séance ElementR, fondée sur le package R [spatstat.linnet](https://cran.r-project.org/web/packages/spatstat.linnet/index.html)
 de la famille [spatstat](https://cran.r-project.org/web/packages/spatstat/index.html).