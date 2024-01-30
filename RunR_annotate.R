#R analysis bfromm

library(seqinr)


getPairs<-function(file_in, seqMin,seqMax){

L1<-file_in[,3]-file_in[,2]
L2<-file_in[,9]-file_in[,8]
Sel<-which(L1>seqMin&L1<seqMax&L2>seqMin&L2<seqMax)

Pair<-file_in[Sel,13]

OutHist<-table(factor(Pair,levels=1:20))
return(OutHist)
}

formals(getPairs)[[2]]<-25
formals(getPairs)[[3]]<-32

siPairs<-function(file_in, seqSel){

L1<-file_in[,3]-file_in[,2]
L2<-file_in[,9]-file_in[,8]

Sel<-which(L1==seqSel&L2==seqSel)
Over3<-abs(file_in[Sel,3]-file_in[Sel,5])
OutHist<-table(factor(Over3, levels=1:10))
return(OutHist)
}

pairfiles<-list.files(pattern="pair_")

PairList<-list()
for(i in 1:length(pairfiles)){
PairList[[i]]<-read.table(pairfiles[i], sep="\t",stringsAsFactors=F)
OverallDist<-PairList[[i]][,3]-PairList[[i]][,2]
TabLength<-table(factor(OverallDist,levels=18:32))

barplot(TabLength/sum(TabLength), ylim=c(0,0.3))
#dev.copy(pdf, paste(pairfiles[i],"_overalllength.pdf",sep=""))
#dev.off()
piPairs<-getPairs(PairList[[i]])
barplot(piPairs/sum(piPairs),ylim=c(0,0.7))
#dev.copy(pdf, paste(pairfiles[i], "_overlapPi.pdf",sep=""))
#dev.off()

Pairsi<-siPairs(PairList[[i]],21)
barplot(Pairsi/sum(Pairsi), ylim=c(0,0.5))
#dev.copy(pdf, paste(pairfiles[i], "_3overhangsi.pdf",sep=""))
#dev.off()


}



#we can also do this specifically for siRNAs that intersect with protein coding genes. 

annotatedFiles<-list.files(pattern="annot")
AnnFileList<-vector("list", length=length(annotatedFiles))

annotatedFiles<-annotatedFiles[-1]
for(i in 1:length(annotatedFiles)){
AnnFileList[[i]]<-read.table(annotatedFiles[i], sep="\t",stringsAsFactors=F)
WFile<-AnnFileList[[i]]
SenseOverlap<-which((WFile[,6]=="+"&WFile[,20]=="+")|(WFile[,6]=="-"&WFile[,20]=="-"))
AntOverlap<-which((WFile[,6]=="+"&WFile[,20]=="-")|(WFile[,6]=="-"&WFile[,20]=="+"))

SenseTab<-table(factor(WFile[SenseOverlap,3]-WFile[SenseOverlap,2],levels=18:32))
AntTab<-table(factor(WFile[AntOverlap,3]-WFile[AntOverlap,2],levels=18:32))
barplot(SenseTab/(sum(SenseTab)+sum(AntTab)), ylim=c(0,0.3))
dev.copy(pdf, paste(annotatedFiles[i], "_lengthSense.pdf",sep=""))
dev.off()

barplot(AntTab/(sum(AntTab)+sum(SenseTab)),ylim=c(0,0.3))
dev.copy(pdf, paste(annotatedFiles[i], "_lengthAnt.pdf",sep=""))
dev.off()

AllOv<-c(SenseOverlap,AntOverlap)
piPairs<-getPairs(WFile[AllOv,])
barplot(piPairs/sum(piPairs),ylim=c(0,0.7))
dev.copy(pdf, paste(annotatedFiles[i], "_piOverlap.pdf",sep=""))
dev.off()

}

#how clustered?

genomeLengths<-list()

pairgenomes<-list.files(pattern="_genom")
pairgenomes<-pairgenomes[c(1,2,2,2,4,4,5,6)]
for(i in 1:length(pairgenomes)){
Gn<-read.fasta(pairgenomes[i])

genomeLengths[[i]]<-sapply(Gn,length)

}

clusterPi<-function(data_in,genomeLengths,minL,maxL){
piSel<-data_in[which(data_in[,3]-data_in[,2]>minL&data_in[,3]-data_in[,2]<maxL),]
piTab<-table(piSel[,1])
q<-match(names(piTab),names(genomeLengths))
piTabN<-piTab/genomeLengths[q]
piTab<-piTab[order(piTabN, decreasing=T)]
piTabN<-piTabN[order(piTabN,decreasing=T)]

return(cbind(piTab,piTabN))
}

AllClust<-list()
for(i in 1:length(PairList)){
AllClust[[i]]<-clusterPi(PairList[[i]],genomeLengths[[i]],25,32)



}


plot(seq(0,1, length=nrow(AllClust[[1]])),cumsum(AllClust[[1]][,1])/sum(AllClust[[1]][,1]),type="l", lwd=2, col="black",ylab="fract_piRNA",xlab="fract_contig")


ColIn<-c("black","red","red","red","blue","blue","grey","brown")
for(i in 2:8){
points(seq(0,1,length=nrow(AllClust[[i]])),cumsum(AllClust[[i]][,1]/sum(AllClust[[i]][,1])),type="l", lwd=2, col=ColIn[i])




}

dev.copy(pdf, "piDensClust_twoNorm.pdf")
dev.off()

plot(cumsum(AllClust[[1]][,1])/sum(AllClust[[1]][,1]),type="l", lwd=2, col="black",ylab="fract_piRNA",xlab="N_contig")


ColIn<-c("black","red","red","red","blue","blue","grey","brown")
for(i in 2:8){
points(cumsum(AllClust[[i]][,1]/sum(AllClust[[i]][,1])),type="l", lwd=2, col=ColIn[i])




}
dev.copy(pdf, "piDensClust_oneNorm.pdf")
dev.off()


plot(log2(cumsum(AllClust[[1]][,1])),type="l", lwd=2, col="black",ylab="fract_piRNA",xlab="N_contig",ylim=c(10,23))


ColIn<-c("black","red","red","red","blue","blue","grey","brown")
for(i in 2:8){
points(log2(cumsum(AllClust[[i]][,1])),type="l", lwd=2, col=ColIn[i])




}
dev.copy(pdf, "piDensClust_NoNorm.pdf")
dev.off()
#density analysis Rdata file.


clusterWin<-function(data_in,minL,maxL){
data_in<-data_in[which(data_in[,ncol(data_in)]>0),]
piSel<-data_in[which(data_in[,6]-data_in[,5]>minL&data_in[,6]-data_in[,5]<maxL),]
piNam<-paste(piSel[,1],piSel[,2],piSel[,3], sep=":")
piTab<-table(piNam)

piTab<-piTab[order(piTab, decreasing=T)]


return(piTab)
}


#compare real clustering to the random clustering

RandomFiles<-list.files(pattern="RandomOverlap")

windowFiles<-list.files(pattern="window25")
Real90All<-c()
for(i in 1:length(RandomFiles)){

RndTemp<-read.table(RandomFiles[i], sep="\t", stringsAsFactors=F)

RndClust<-clusterWin(RndTemp, 25,32)
WinTemp<-read.table(windowFiles[i], sep="\t", stringsAsFactors=F)

WinClust<-clusterWin(WinTemp,25,32)
XaxNorm1<-seq(0,1,length=length(RndClust))
XaxNorm2<-seq(0,1,length=length(WinClust))
RndCumSum<-cumsum(RndClust)/max(cumsum(RndClust))
Rnd90<-which(RndCumSum>0.9)[1]


WinCumSum<-cumsum(WinClust)/max(cumsum(WinClust))
Real90<-which(WinCumSum>0.9)[1]
RndEq<-RndCumSum[Real90]
Real90All<-rbind(Real90All, c(Rnd90,Real90,RndEq))
#plot(XaxNorm1,RndCumSum,ylab="fraction_piRNAs",xlab="fraction_windows",type="l", lwd=3,col="grey")
#points(XaxNorm2,WinCumSum, col="black", type="l", lwd=2)
#dev.copy(pdf, paste("realvNorm_",windowFiles[i], ".pdf", sep=""))
#dev.off()


}


barplot(Real90All)







