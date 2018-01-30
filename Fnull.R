


### Functional null model


spxplot<-rbind(S1=c(A=1,B=1,C=1,D=0,E=0),
      S2=c(A=1,B=1,C=0,D=1,E=0),
      S3=c(A=1,B=0,C=0,D=1,E=1))

spxtrait<-rbind(A=c(t1=1,t2=1),
                B=c(t1=2,t2=2),
                C=c(t1=3,t2=2),
                D=c(t1=3,t2=3),
                E=c(t1=2,t2=1))

coords<-c(1,2,3)


### 

library(vegan)
library(ecodist)
library(betapart)
library(FD)

jac<-vegdist(spxplot,"jac")
geodist<-dist(coords)

plot(jac~geodist)

functional.beta.pair(as.matrix(spxplot),(spxtrait),"jaccard")


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


func.pair<-functional.beta.pair(x=comm.test, traits=traits.test, index.family = "jaccard" )


taxo.pair<-beta.pair(x=comm.test,index.family = "jaccard" )

cor(taxo.pair$beta.jtu,func.pair$funct.beta.jtu)

ver<-{}
for(i in 1:1000){
  tryCatch({
    traits.test.random<-apply(traits.test,2,function(x)sample(x))
    func.pair.random<-functional.beta.pair(x=comm.test, traits=traits.test.random, index.family = "jaccard" )
    
    ver[i]<-cor(taxo.pair$beta.jac,func.pair.random$funct.beta.jac)
    
      }, error=function(e){})
}

hist(ver)




#######

## With real data


otto3<-read.table("otto3.txt",header=TRUE)
traits<-read.table("functional traits.txt",header=TRUE)

sum(!row.names(traits)%in%colnames(otto3))# must be zero
sum(!colnames(otto3)%in%row.names(traits))# must be zero


##

gower.dist<-gowdis(traits)
pcoa<-cmdscale(gower.dist,5)

func.pair<-functional.beta.pair(x=otto3, traits=pcoa[,1:2], index.family = "jaccard" )

taxo.pair<-beta.pair(x=comm.test,index.family = "jaccard" )

traits.test.random<-apply(traits.test,2,function(x)sample(x



cor(taxo.pair$beta.jtu,func.pair$funct.beta.jtu)

ver<-{}
for(i in 1:1000){
  tryCatch({
    traits.test.random<-apply(traits.test,2,function(x)sample(x))
    func.pair.random<-functional.beta.pair(x=comm.test, traits=traits.test.random, index.family = "jaccard" )
    
    ver[i]<-cor(taxo.pair$beta.jac,func.pair.random$funct.beta.jac)
    
  }, error=function(e){})
}

hist(ver)













  