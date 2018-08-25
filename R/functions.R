# Functions definitions

#' Identificar la presencia de ensambles a partir de análisis de agrupamiento (cluster)
#'
#' @param epifauna abundancia de especies
#' 
#' @library vegan, gclus, cluster, RColorBrewer,labdsv
#'
#' @return análisis cluster de agrupamiento para identificar la presencia de ensambles
#' @export
#'
#'Data accessibility
#'
#' @examples
#' 
Cluster <- function(epifauna) 
{
  epifauna <- data.frame(epifauna[,2:14], row.names = epifauna$Station)
  epifauna.hell <- decostand(epifauna, "hellinger")
  epifauna.hell.bc <- vegdist(epifauna.hell, method = "bray")
  clus2 <- hclust(epifauna.hell.bc, method = "ward.D2")
  den2 <- as.dendrogram(clus2)
  H <- plot(den2, xlab = "Stations", ylab = "Bray-Curtis dissimilarity", nodePar = list(pch = c(1, NA), lab.cex = 0.7))
  return(H)
}


# epifauna <- data.frame(epifauna[,2:14], row.names = epifauna$Station) 
# 'Station' como aparte de la tabla
# Chequear si 2:14 se podría reemplazar por algo tipo i:n para decirle que tome desde la segunda hasta el final de los datos


# epifauna.hell <- decostand(epifauna, "hellinger")
# Transformación de Hellinger

# epifauna.hell.bc <- vegdist(epifauna.hell, method = "bray")
# Bray Curtis

# clus2 <- hclust(epifauna.hell.bc, method = "ward.D2")
# den2 <- as.dendrogram(clus2)
# plot(den2, xlab = "Stations", ylab = "Bray-Curtis dissimilarity", nodePar = list(pch = c(1, NA), lab.cex = 0.7))
# return(plot)
# Cluster "Ward's"