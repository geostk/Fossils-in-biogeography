#########################dependencies#####################################
require(ape)
require(phangorn)
##########################################################################
############################Iputs#########################################
wd = "C:\\Users\\Ruud\\Desktop\\R working Directory\\adding fossils to tree"
setwd(wd)
trees_sub <- read.nexus("100trees_23041989.nex") # load in tree dataset
fossil_dir = "Fossil_lists\\" #set directory where the fossil files are
##########################################################################
##########################################################################

place.fossil<-function(file,tr,fossil_nr){ 
  
  #load in dependencies
  require(ape)
  require(phangorn)
  
  
  #load in fossil data
  fossil <-list()
  fossilsource = file
  suppressWarnings(fossil$general<-read.table(fossilsource,nrows=1, header = TRUE,stringsAsFactors = FALSE))
  suppressWarnings(fossil$decendants<-scan(fossilsource,skip=2,what="character"))
  
  print(paste(tree_nr," ",fossil_nr, sep = ""))
  
  #Generate test tree
  #trstr = "((orang:5,(numnum:3,bumbum:3):2):6,(gorilla:7,(chimp:6,human:6):1):4);"
  #tr = read.tree(text=trstr)
  
  #set min/max fossilage
  minage = as.numeric(fossil$general[,2])
  maxage = as.numeric(fossil$general[,3])
  
  #get list of species (decendants of the fossil)
  list <- fossil$decendants
  fossilage = runif(1, minage, maxage) # select random fossilage between the minimum and the naximumage
  
  #print(paste(tree_nr," ",fossil_nr," B", sep = ""))
  
  #get all edges that belong to the clade 
  #if there is only one tip in the list:
  if (length(list)==1){
    tipnr<-grep(list,tr$tip.label) #find the number of the singular tip
    edgenode <- which(tr$edge[,2] %in% tipnr) # select the edge connected to the tip
  } else{ #when there are multiple tips
    list_numbers<-which(tr$tip.label %in% list)
    MRCA<-mrca.phylo(tr,list_numbers)
    #MRCA<-getMRCA(tr,list) #get MRCA
    decendant_nodes <- Descendants(tr,MRCA,"all") #find the decendant nodes
    clade_tips <- decendant_nodes[decendant_nodes <= length(tr$tip.label)] #find the names of the decendant tips
    edgenode <-which.edge(tr,clade_tips) #get all edges that connect the tips
    #for a stem fossile you only the sten
    if (fossil$general[4] == "stem"){
      edgenode <-which(tr$edge[,2] %in% MRCA)
    }
    #in the case of both its all crown edges + stem
    if (fossil$general[4] == "both"){
      edgenode[length(edgenode)+1] <-which(tr$edge[,2] %in% MRCA)
    }
  }
  
  #print(paste(tree_nr," ",fossil_nr," C", sep = ""))
  
  #find all branches ages
  branchtimes = branching.times(tr)
  branchages <- tr$edge
  for (y in 1:length(tr$edge)){
    if (tr$edge[y] > length(tr$tip.label)){
      branchages[y] = branchtimes[[tr$edge[y]-length(tr$tip.label)]] 
    }else{
      branchages[y] = 0
    }
  }
  
  #make sure the fossil age isnt older than  the node age
  
  if (max(branchages[edgenode,])<fossilage){
    if (max(branchages[edgenode,]) >= minage){
      fossilage = runif(1, minage, max(branchages[edgenode,]))
    } else{
        fossilage = max(branchages[edgenode,])-0.001
        warning("Fossil min age is older than the maximum nodeage of the clade.
Fossilage is set to equal the clade's maxage.")
    }
  }
  
  #print(paste(tree_nr," ",fossil_nr," D", sep = ""))
  
  # put starting and finishing ages of selected branches in "between" object
  goodnodes <- vector()
  o = 1
  between <-matrix(branchages[edgenode,],ncol = 2)
  #Calculate the age at which the fossil connects
  if (fossil$general[4] != "stem"){
    fossilconnectionage <- runif(1, fossilage, max(between))
  } else{
    if(min(between)>fossilage){
      fossilconnectionage <- runif(1, min(between), max(between))
    } else{
      fossilconnectionage <- runif(1, fossilage, max(between))
    }
  }
  #check if the fossils connection age is between an edge's min/max age
  for (x in 1:length(between[,1])){
    if (between[x,1] > fossilconnectionage & between[x,2] < fossilconnectionage){
      goodnodes[o] <- edgenode[x]
      o = o+1
    }
  }
  
  #print(paste(tree_nr," ",fossil_nr," E", sep = ""))
  
  #randomly select the branch to place the fossil
  randombranch<-goodnodes[sample(1:length(goodnodes),1)]
  #select upper node (weird name for object I know)
  tip1 = tr$edge[randombranch,1]
  #select lower node (weird name for object I know)
  tip2 = tr$edge[randombranch,2]
  
  #############################################################################################
  #Adding fossils to tree#######################################################################
  #############################################################################################
  #############################################################################################
  
  
  ##################################Some stuff i use for tests#####################
  #add_fossil <- function(tr=tr,tip1=tip1,tip2=tip2,fossilage=fossilage,fossilconnectionage=fossilconnectionage){
  #library(ape)
  
  #trstr = "((orang:5,(numnum:3,bumbum:3):2):6,(gorilla:7,(chimp:6,human:6):1):4);"
  #tr = read.tree(text=trstr)
  
  
  #tip1 = 8
  #tip2 =9
  #fossilage = 2
  #fossilconnectionage = 4
  ##################################################################
  ###################################################################
  
  tr2 = tr
  
  #calculate branchlengths of the new branches
  branchtimes = branching.times(tr)
  fossilconnection_top_branchlength = branchtimes[[tip1-length(tr$tip.label)]]-fossilconnectionage
  
  #print(paste(tree_nr," ",fossil_nr," F", sep = ""))
  
  if(tip2 > length(tr$tip.label)){
    fossilconnection_lower_branchlength = fossilconnectionage-branchtimes[[tip2-length(tr$tip.label)]]
  }
  fossil_branchlength = fossilconnectionage-fossilage
  
  
  #expand the matrix of edges for the new tree
  tr2$edge = matrix(ncol = 2,nrow = length(tr$edge[,1])+2)
  tr2$tip.label = append(tr$tip.label,fossil$general[[1]],after=0)
  tr2$Nnode = tr$Nnode+1
  
  #fix edge numbers
  tempedge = tr$edge+1
  #find the branch on which a the fossil is added
  change1 = which(tr$edge[,2] %in% tip2)
  
  #print(paste(tree_nr," ",fossil_nr," G", sep = ""))
  
  #add the fossil to the tree if it is added on a branch that leads to a tip
  if(tip2 <= length(tr$tip.label)){
    #print(paste(tree_nr," ",fossil_nr," Line 1", sep = ""))
    for (i in 1:length(tr2$edge[,2])){
      #before the branch keep old numbering +1
      if (i < change1[1]){
        tr2$edge[i,]=tempedge[i,]
        tr2$edge.length[i]=tr$edge.length[i]
      }
      if (i == change1[1]){
        tr2$edge[i,]=c(tip1+1,tip1+2)
        tr2$edge.length[i]=fossilconnection_top_branchlength
      }
      #add the edges
      if (i == (change1[1]+1)){
        tr2$edge[i,] = c((tip1+2),1) 
        tr2$edge.length[i]=fossil_branchlength
      }
      if (i == (change1[1]+2)){
        tr2$edge[i,] = c(tip1+2,tip2+1) 
        tr2$edge.length[i]=fossilconnectionage
      }
      
      
      #print(paste(tree_nr," ",fossil_nr," H", sep = ""))
      
      #after contuniue old numbering +2 (with the exeption of tips)
      if (i > (change1[1]+2)){ # after the fossil inclusion
        if (tempedge[(i-2),2] > length(tr2$tip.label) & tempedge[(i-2),1] > tip1+1){
          tr2$edge[i,] = tempedge[(i-2),]+1 
          tr2$edge.length[i]=tr$edge.length[i-2]
        } 
        if (tempedge[(i-2),2] <= length(tr2$tip.label)){ # don't add an extra node number to lower tip when dealing with tips
          tr2$edge[i,] = c(tempedge[(i-2),1]+1,tempedge[(i-2),2]) 
          tr2$edge.length[i]=tr$edge.length[i-2]
        }
        if (tempedge[(i-2),1] <= tip1+1){ #if the parent node is a smaller number than tip1, we don't need to add the extra +1
          if (tempedge[(i-2),2] <= length(tr$tip.label)+1){ # don't add an extra node number to lower tip when dealing with tips
            tr2$edge[i,] = c(tempedge[(i-2),1],tempedge[(i-2),2])
          } else {
            tr2$edge[i,] = c(tempedge[(i-2),1],tempedge[(i-2),2]+1)
          }
          tr2$edge.length[i]=tr$edge.length[i-2]
        }
      }
    }
  } else { #add the fossil to the tree if it is added on an internal branch
    #print(paste(tree_nr," ",fossil_nr," Line 2", sep = ""))
    for (i in 1:length(tr2$edge[,2])){
      #before the branch keep old numbering +1
      if (i < change1[1]){
        tr2$edge[i,]=tempedge[i,]
        tr2$edge.length[i]=tr$edge.length[i]
      }
      if (i == change1[1]){
        tr2$edge[i,]=tempedge[i,]
        tr2$edge.length[i]=fossilconnection_top_branchlength
      }
      #add the edges
      if (i == (change1[1]+1)){
        tr2$edge[i,] = c((tip2+1),1) 
        tr2$edge.length[i]=fossil_branchlength
      }
      if (i == (change1[1]+2)){
        tr2$edge[i,] = c(tip2+1,tip2+2) 
        tr2$edge.length[i]=fossilconnection_lower_branchlength
      }
    
      #after contuniue old numbering +2 (with the exeption of tips)
      if (i > (change1[1]+2)){
        if (tempedge[(i-2),2] > length(tr2$tip.label)){
          tr2$edge[i,] = tempedge[(i-2),]+1 
          tr2$edge.length[i]=tr$edge.length[i-2]
        } 
        if (tempedge[(i-2),2] <= length(tr2$tip.label)){
          tr2$edge[i,] = c(tempedge[(i-2),1]+1,tempedge[(i-2),2]) 
          tr2$edge.length[i]=tr$edge.length[i-2]
        }
        if (tempedge[(i-2),1] <= tip1+1){
          if (tempedge[(i-2),2] <= length(tr$tip.label)+1){
            tr2$edge[i,] = c(tempedge[(i-2),1],tempedge[(i-2),2])
          } else {tr2$edge[i,] = c(tempedge[(i-2),1],tempedge[(i-2),2]+1)
          }
          tr2$edge.length[i]=tr$edge.length[i-2]
        }
      }
    }
  }
  #print(paste(tree_nr," ",fossil_nr," I", sep = ""))
  return(tr2)
  if (length(tr2$edge.length[tr2$edge.length<0])>0){
    print(paste(tree_nr," ",fossil_nr," NEGATIVE_BRANCHES", sep = ""))
  }
}

#function that places all fossils within the fossil directory in the tree
place.multiple.fossils <-function(fossil_dir,tr){
  #tr <- trees_sub[[11]] 
  #tr <- ladderize(tr,FALSE)
  fossillist<-dir(fossil_dir)
  write.nexus(file="temp.nex",tr) #writes a temporary file that it reloads later, did this because I got some weitd errors regarding the numbering. However I'm not sure it's neccisairy anymore
  for (i in 1:length(fossillist)){ # loop over all fossil files
    fossil_nr <- paste("fossil ",i,sep="")
    #print(fossil_nr)
    tr <- read.nexus("temp.nex")
    #print(paste(tree_nr," ",fossil_nr," O2", sep = ""))
    fossilsource <- paste(fossil_dir,fossillist[i],sep="") # get the exact name of the fossil file (from all files in that directory)
    #print(paste(tree_nr," ",fossil_nr," O3", sep = ""))
    tr<-place.fossil(file = fossilsource,tr=tr,fossil_nr=fossil_nr) # the previous finctions that places a fossil
    #print(paste(tree_nr," ",fossil_nr," O4", sep = ""))
    #plot(tr)
    #print(paste(tree_nr," ",fossil_nr," O5", sep = ""))
    write.nexus(file="temp.nex",tr) # replace temp.nex with the new tree
    #print(paste(tree_nr," ",fossil_nr," O6", sep = ""))
    #i=i+1
    #print(i-1)
  }
  
  return(tr)
  #print(paste(tree_nr," ",fossil_nr," return tree", sep = ""))
}


trees <- rmtree(length(trees_sub), 2) # create random multiphylo file to fill witht he new tree later
for (i in 1:length(trees_sub)){ #for each of the trees do the following:
  tree_nr <- paste("tree_number_",i,sep="")
  #print(tree_nr)
  tr <- ladderize(trees_sub[[i]],FALSE) #ladderize the tree
  tr<-place.multiple.fossils(fossil_dir=fossil_dir,tr=tr) #place all fossils in each tree and add it to multiphylo object
  trees[[i]]<-tr #place tree in the multiphylo object (was the original output, but kept it there because it is ocasionaly handy for bugtracking)
  tree_name = paste(tree_nr,".tree",sep="")
  if(length(tr$edge.length[tr$edge.length<0])>0){
    tree_name = paste(tree_nr,"_NEGATIVE_BRANCHES",".tree",sep="")
  }
  write.tree(tr, file = paste(tree_nr,".tree",sep=""))
  rm(tree_nr,tr,tree_name) 
  gc(verbose = FALSE)
  
  #pdf(file = tree_name, width= 100, height = 100)
  #plot(tr)
  #nodelabels()
  #edgelabels()
  #tiplabels()
  #axisPhylo()
  #dev.off()
}


