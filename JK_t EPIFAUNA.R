# SPATIAL ANALYSIS OF BENTHIC BIODIVERSITY, SAN JORGE GULF ARGENTINA
# Julieta Kaminsky
# 01/2017

# R reference
citation()
citation("vegan")

# packages
library(vegan)
library(gclus)
library(cluster)
library(RColorBrewer)
library(labdsv)
library(readxl)

# Taxonomic diversity

# Import .xls
library(readxl)
t.2000 <- read_excel(file.choose(),1)
head(t.2000)
str(t.2000)
t.2000 <- data.frame(t.2000[,2:52], row.names = t.2000$Station)

# Transformation Hellinger
library(vegan)
t.2000.hell <- decostand(t.2000, "hellinger")

# Transformation fourth-root
#t.2000.trans <- sqrt(sqrt(t.2000))

# Bray-Curtis distance
library(vegan)
t.2000.hell.bc <- vegdist(t.2000.hell, method = "bray")
t.2000.hell.bc

#par(mar = (c(3,4,1,0.3)))
#layout(matrix(c(1:2), 2, 1))

# Cluster "Ward's"
clus2 <- hclust(t.2000.hell.bc, method = "ward.D2")
den2 <- as.dendrogram(clus2)
plot(den2, xlab = "Stations", ylab = "Bray-Curtis dissimilarity", nodePar = list(pch = c(1, NA), lab.cex = 0.7))
rect.hclust(clus2, h = 0.8 , border = 2:6)

# nMDS
(t.2000.mds2 <- metaMDS(t.2000.hell.bc, k = 2))
(t.2000.mds3 <- metaMDS(t.2000.hell.bc, k = 3))
(t.2000.mds4 <- metaMDS(t.2000.hell.bc, k = 4))
# With 2 Stress is 0.1631

#layout(matrix(c(1:3), 1, 3))
#stressplot(t.14.mds2, t.14.hell.bc, main = "2 dimensions")
#stressplot(t.14.mds3, t.14.hell.bc, main = "3 dimensions")
#stressplot(t.14.mds4, t.14.hell.bc, main = "4 dimensions")
#layout(1)

# Fig in 2 dimensions
#plot(t.2000.mds2, choices = c(1, 2), display = "sites", type = "text")
# Stress 0.185

# nMDS + Cluster
# Add colours from a Cluster result to an nMDS plot
clus2
clus2.groups <- cutree(clus2, h = 0.8)
grp.lev <- levels(factor(clus2.groups))

# Combination with nMDS result
sit.sc <- scores(t.2000.mds2)
p <- ordiplot(sit.sc, type="n")
for (i in 1:length(grp.lev)) {
  points(sit.sc[clus2.groups==i,], pch=(14+i), cex=2, col=i+1)
}
text(sit.sc, row.names(t.2000), pos=4, cex=0.7)

# ordicluster(p, clus2, col="drak grey")
legend(locator(1), c("Group A", "Group ", "Group ", "Group ", "Group ", "Group "), pch=14+c(1:length(grp.lev)), col=1+c(1:length(grp.lev)), pt.cex=1.5, cex = 1)
text(locator(1), label="Stress = 0.16", col="black", cex=0.8)

# SIMPER
library(readxl)
t.2000.simper <- read_excel(file.choose(), 1)
head(t.2000.simper)
str(t.2000.simper)

groupes <- t.2000.simper$Group
t.2000.simper <- data.frame(t.2000.simper[,3:53], row.names = t.2000.simper$Station)

library(vegan)
sim <- simper(t.2000.simper, groupes, permutations = 9999)
summary(sim)
sim

# Diversity indices 
# Alpha diversity by station
library(vegan)
t.2000.alpha <- specnumber(t.2000)
t.2000.alpha <- data.frame(t.2000.alpha)

# Gamma regional
col_03 <- colSums(t.2000)
library(vegan)
t.2000.gamma <- specnumber(col_03)
t.2000.gamma <- data.frame(t.2000.gamma)


# Alpha diversity by group


# Gamma group

# Group A
library(readxl)
t.2000.A <- read_excel(file.choose(), 1)
t.2000.A <- data.frame(t.2000.A[,3:53], row.names = t.2000.A$Station)

col_2000.A <- colSums(t.2000.A)
library(vegan)
t.2000.A.gamma <- specnumber(col_2000.A)
t.2000.A.gamma <- data.frame(t.2000.A.gamma)

# Group B
library(readxl)
t.2000.B <- read_excel(file.choose(), 1)
t.2000.B <- data.frame(t.2000.B[,3:53], row.names = t.2000.B$Station)

col_2000.B <- colSums(t.2000.B)
library(vegan)
t.2000.B.gamma <- specnumber(col_2000.B)
t.2000.B.gamma <- data.frame(t.2000.B.gamma)


# Group C
library(readxl)
t.2000.C <- read_excel(file.choose(), 1)
t.2000.C <- data.frame(t.2000.C[,3:53], row.names = t.2000.C$Station)

col_2000.C <- colSums(t.2000.C)
library(vegan)
t.2000.C.gamma <- specnumber(col_2000.C)
t.2000.C.gamma <- data.frame(t.2000.C.gamma)


# Group D
library(readxl)
t.2000.D <- read_excel(file.choose(), 1)
t.2000.D <- data.frame(t.2000.D[,3:53], row.names = t.2000.D$Station)

col_2000.D <- colSums(t.2000.D)
library(vegan)
t.2000.D.gamma <- specnumber(col_2000.D)
t.2000.D.gamma <- data.frame(t.2000.D.gamma)

# Group E
library(readxl)
t.2000.E <- read_excel(file.choose(), 1)
t.2000.E <- data.frame(t.2000.E[,3:53], row.names = t.2000.E$Station)

col_2000.E <- colSums(t.2000.E)
library(vegan)
t.2000.E.gamma <- specnumber(col_2000.E)
t.2000.E.gamma <- data.frame(t.2000.E.gamma)

# Group F
library(readxl)
t.2000.F <- read_excel(file.choose(), 1)
t.2000.F <- data.frame(t.2000.F[,3:53], row.names = t.2000.F$Station)

col_2000.F <- colSums(t.2000.F)
library(vegan)
t.2000.F.gamma <- specnumber(col_2000.F)
t.2000.F.gamma <- data.frame(t.2000.F.gamma)

# Beta group