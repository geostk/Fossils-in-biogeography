library(ape)
library(geiger)

geog <- read.table(file = "gent_wFossils_04052015.geog", row.names = 1,skip = 2,)

for (i in 1:100){
  tr <- read.tree(paste("tree_number_",i,".tree",sep=""))
  cleaned <- treedata(tr,geog)
  tr <- cleaned$phy
  write.tree(tr, file = paste("tree_number_",i,".tree",sep=""))
}

