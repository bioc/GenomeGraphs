##################################
#Scripts to produce dummy data
#
##################################

library(GenomeGraphs)

mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")

probestart = seq(180292097,180492096, by = 1000)
p1 = seq(180292097,180492096, by = 334)
p2 = seq(180292097,180492096, by = 273)
probestart = c(probestart, p1, p2)
ord = order(probestart)
probestart = probestart[ord]

cn = rnorm(sum(probestart < 180320000) ,mean = 1, sd = 0.3)
cn1 = rnorm((sum(probestart >=  180320000) - sum(probestart >=  180450000)) ,mean = 2.8, sd = 0.3)
cn2 = rnorm(sum(probestart > 180450000) ,mean = 1, sd = 0.3)
cn = as.matrix(c(cn, cn1,cn2))

minbase = 180292097 
maxbase = 180492096


segments = list(c(1,2.8,1))
segStart = list(c(180292097, 180320000,180450000))
segEnd = list(c(180320000,180450000, 180492096))

genesplus  =new("GeneRegion", start = minbase, end = maxbase, strand = "+", chromosome = "3", biomart=mart)
genesmin  =new("GeneRegion", start = minbase, end = maxbase, strand = "-", chromosome = "3", biomart=mart)

exo=NULL
for(i in 1:length(genesplus@ens[,1])){
 exo = c(exo,sample(genesplus@ens[i,4]:genesplus@ens[i,5],4))
}
ord = order(exo)
exonProbePos = exo[ord]

int1 = rnorm(sum(exo < 180320000), mean = 3, sd = 1) 
int2 = rnorm((sum(exo >=  180320000) - sum(exo >=  180450000)) ,mean = 7, sd = 0.5)
int3 = rnorm(sum(exo > 180450000) ,mean = 2.7, sd = .8)

intensity = as.matrix(c(int1,int2,int3))

save(intensity, exonProbePos, cn, probestart, segments, segStart, segEnd, yeastCons1, file = "dummyData.rda")

####################

