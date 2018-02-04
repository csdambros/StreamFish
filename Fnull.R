##### Functional diversity of stream fishes in Brazil - Analysis
##### Created by CSDambros in Feb 1st 2018


### Load packages

library(vegan)
library(ecodist)
library(betapart)
library(FD)

### Source functions created or modified
source("newbeta.part.core.R")

###
traits.test<-cbind( c(1,1,1,2,2,3,3,4,4,5,5) , c(1,2,4,1,2,3,5,1,4,3,5) )
dimnames(traits.test)<-list(paste("sp",1:11,sep="") , c("Trait 1","Trait 2") ) 

comm.test<-matrix(0,4,11,dimnames=list( c("A","B","C","D") , paste("sp",1:11,sep="") ) )
comm.test["A",c(1,2,4,5)]<-1
comm.test["B",c(1,3,8,9)]<-1
comm.test["C",c(6,7,10,11)]<-1
comm.test["D",c(2,4,7,9)]<-1


plot(5,5,xlim=c(0,6), ylim=c(0,6), type="n", xlab="Trait 1",ylab="Trait 2")
points(traits.test[,1],traits.test[,2], pch=21,cex=1.5,bg="black")
rect(1,1,4,4, col="#458B0050", border="#458B00") ; text(2.5,2.5,"B",col="#458B00",cex=1.5)	
polygon(c(2,1,3,4), c(1,2,5,4), col="#DA70D650", border="#DA70D6") ; 
text(2.5,3,"D",col="#DA70D6",cex=1.5)	
rect(1,1,2,2, col="#FF000050" , border="#FF0000") ; text(1.5,1.5,"A",col="#FF0000",cex=1.5)	
rect(3,3,5,5, col="#1E90FF50", border="#1E90FF") ; text(4,4.2,"C",col="#1E90FF",cex=1.5)	


dir.create("verts")
dir.create("resu")

func.pair<-functional.beta.pair(x=comm.test, traits=traits.test, index.family = "jaccard",prefix = "verts/obs")

taxo.pair<-beta.pair(x=comm.test,index.family = "jaccard" )

# Randomize traits and calculate functional similarity
func.pair.random<-lapply(1:100,functional.beta.pair.random,comm=comm.test,traits=traits.test,index.family="jaccard",prefix="verts/rtest1.",gower=FALSE)

func.pair2.random<-lapply(1:100,functional.beta.pair2.random,comm=comm.test,traits=traits.test,index.family="jaccard",prefix="verts/rtest1.",gower=FALSE)

func.pair2.random<-lapply(func.pair2.random,function(x){lapply(x,function(x){as.dist(t(x))})})


length(func.pair.random)
length(func.pair2.random)

# Remove matrices where error occurred (null)
func.pair.random<-func.pair.random[unlist(lapply(func.pair.random,function(x)!is.null(x)))]
func.pair2.random<-func.pair2.random[unlist(lapply(func.pair2.random,function(x)length(x)>0))]

# Check how many aleatorizations remained
length(func.pair.random)
length(func.pair2.random)

# Correlate taxonomic and functional similarity expected under null model
hist(unlist(lapply(func.pair.random, function(x)cor(x$funct.beta.jac,taxo.pair$beta.jac))))
hist(unlist(lapply(func.pair2.random, function(x)cor(x$jac,taxo.pair$beta.jac))))

#### Test in parallel

library(snow)

criscl<-makeCluster(20)

clusterApply(criscl,1:15,function(x){Sys.sleep(2);system("ifconfig | grep '192.168'")})

clusterExport(criscl,list("functional.beta.pair.random","functional.betapart.core","functional.beta.pair","traits.test","comm.test"))

func.pair.random<-clusterApply(criscl,1:100,functional.beta.pair.random,comm=comm.test,traits=traits.test,index.family="jaccard",prefix="verts/rtest1.",gower=FALSE)


stopCluster(cl)

length(func.pair.random)

# Remove matrices where error occurred (null)
func.pair.random<-func.pair.random[unlist(lapply(func.pair.random,function(x)!is.null(x)))]

# Check how many aleatorizations remained
length(func.pair.random)

# Correlate taxonomic and functional similarity expected under null model
hist(unlist(lapply(func.pair.random, function(x)cor(x$funct.beta.jac,taxo.pair$beta.jac))))




#######

## With real data


otto3<-read.table("otto3.txt",header=TRUE)
traits<-read.table("functional traits.txt",header=TRUE)

sum(!row.names(traits)%in%colnames(otto3))# must be zero
sum(!colnames(otto3)%in%row.names(traits))# must be zero


## Using with real data

taxo.pair<-beta.pair(x=otto3,index.family = "jaccard" )

dir.create("resu")


gower.dist<-gowdis(traits)
pcoa<-cmdscale(gower.dist,2)

#func.pair<-functional.beta.pair(x=otto3[1:5,], traits=pcoa, index.family = "jaccard",prefix="verts/obs")
func.pair<-functional.beta.pair(x=otto3, traits=pcoa, index.family = "jaccard",prefix="verts/obs")

# Run a single time to check if is working
functional.beta.pair.random(1,comm=otto3[1:2,],traits=pcoa,index.family="jaccard",prefix="rtest1.",gower = FALSE)

functional.beta.pair2.random(1,comm=otto3[1:2,],traits=pcoa,index.family="jaccard",prefix="rtest1.",gower = FALSE)


# Randomize traits and calculate functional similarity
#func.pair.random<-lapply(1:10,functional.beta.pair.random,comm=otto3[1:2,],traits=pcoa,index.family="jaccard",prefix="verts/rtest1.",gower=FALSE)

## Running in parallel

#### Test in parallel

library(snow)

cl<-makeCluster(20)

clusterApply(cl,1:20,function(x){Sys.sleep(2);system("ifconfig | grep '192.168'")})

clusterExport(cl,list("functional.beta.pair.random","functional.beta.pair2.random","functional.betapart.core","functional.beta.pair","traits","otto3","pcoa"))

{

#system.time(
func.pair.random<-clusterApply(cl,1:1000,functional.beta.pair2.random,comm=otto3,traits=pcoa,index.family="jaccard",prefix="verts/rtest1.",gower=FALSE)
#)


#system.time(
#func.pair.random<-clusterApply(cl,21:40,functional.beta.pair.random,otto3[1:10,],traits=pcoa,index.family="jaccard",prefix="verts/rtest1.",gower=FALSE)
#)
  
  
stopCluster(cl)
  
}


#   
# func.pair.random<-clusterApply(cl,1:100,function(i)functional.beta.pair.random(i,comm=otto3[1:15,],traits=traits[,],index.family="jaccard",prefix="verts/rtest1."))

  
  
  
length(func.pair.random)

# Remove matrices where error occurred (null)
func.pair.random<-func.pair.random[unlist(lapply(func.pair.random,function(x)!is.null(x)))]

# Check how many aleatorizations remained
length(func.pair.random)

# Correlate taxonomic and functional similarity expected under null model
hist(unlist(lapply(func.pair.random, function(x)cor(x$funct.beta.jac,taxo.pair$beta.jac))))






############## Ends here



taxo.pair<-beta.pair(x=comm.test,index.family = "jaccard")

cor(taxo.pair$beta.jtu,func.pair$funct.beta.jtu)

func.pair.otto3.random<-list()
ver<-{}


func.pair.random.otto3<-lapply(1:2,functional.beta.pair.random,x=otto3,traits=traits, index.family = "jaccard",prefix = "verts/sim_")



for(i in 1:100){
  tryCatch({
    traits.random<-traits[sample(1:nrow(traits)),]
    rownames(traits.random)<-rownames(traits)
    
    gower.dist<-gowdis(traits.random)
    pcoa<-cmdscale(gower.dist,5)
    
    func.pair<-functional.beta.pair(x=otto3, traits=pcoa[,1:2], index.family = "jaccard",prefix="verts/obs")

    func.pair.random[[i]]<-functional.beta.pair(x=comm.test, traits=traits.test.random, index.family = "jaccard",prefix = paste0("verts/sim_",i))
    
    ver[i]<-cor(taxo.pair$beta.jac,func.pair.random[[i]]$funct.beta.jac)
    
  }, error=function(e){})
}















  