#Example usage: Rscript 0_coor2fasta.R inputfile.vcf or Rscript	0_coor2fasta.R inputfile.bed

require('BSgenome.Hsapiens.UCSC.hg19')
require('rtracklayer')

#process input arguments
args <- commandArgs(trailingOnly = TRUE)


process.simple.offsets.python<-function(DF,prefix="./",window=1000,filterUnique=T,filterCoding=F){
  colnames(DF)[1:6]<-c("A","B","C","D","E","F")
  DF<-DF[!is.na(DF$D),]
  if(filterUnique)
    DF=DF[!duplicated(paste(as.character(DF$C),DF$D,DF$A,DF$B)),]
  DF$D=as.numeric(as.character(as.matrix(DF$D)))
  
  
  halfw=window/2
  DF$C=gsub("chr","",DF$C)
  DF$C=as.character(DF$C)
  DF$C[as.character(DF$C)=="23"]<-"X"
  DF$C[as.character(DF$C)=="24"]<-"Y"
  DF<-DF[as.character(DF$C)!="25",]
  DF<-DF[as.character(DF$C)!="MT",]
  
  #filter by chromosome length
  chrends<-seqlengths(BSgenome.Hsapiens.UCSC.hg19)[match( paste("chr",DF$C,sep=""),names(seqlengths(BSgenome.Hsapiens.UCSC.hg19)))]
  DF=DF[ (as.numeric(as.character(DF$D))-(halfw-1))>0 &  (as.numeric(as.character(DF$D))+halfw)<chrends & (!is.na(chrends)),]
  #filter coding variants
  if(filterCoding){
    require('TxDb.Hsapiens.UCSC.hg19.knownGene')
    library(VariantAnnotation)
    input <- GRanges( seqnames = Rle(paste("chr",DF$C,sep="")),
                      ranges   = IRanges(DF$D, end=DF$D),
                      strand   = Rle(strand("*")) )
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    loc <- locateVariants(input, TxDb.Hsapiens.UCSC.hg19.knownGene, CodingVariants())
    
    DF<-DF[is.na(match(paste(paste("chr",as.character(DF$C),sep=""),DF$D),
                       unique(paste(space(as(loc,"RangedData")),start(loc))))),]
  }
  
  DF.rd1000<-RangedData(ranges=IRanges(start=as.numeric(as.character(DF$D))-(halfw-1),as.numeric(as.character(DF$D))+halfw),space=paste("chr",as.character(DF$C),sep=""),ori=DF$A,mut=DF$B,name=DF$F)

  #retrieve reference genome sequences
  DF.wt1000<-as.character(as.matrix(as.data.frame(getSeq(x=Hsapiens,names=as(DF.rd1000,"GRanges")))))
  
  #write to fasta format
  temp<-DNAStringSet(DF.wt1000)
  names(temp)<-paste(DF.rd1000$ori,DF.rd1000$mut,space(DF.rd1000),start(DF.rd1000),end(DF.rd1000),DF.rd1000$name,sep = "_")
  writeXStringSet(temp,filepath=paste(prefix,".wt",window,".fasta",sep=""))
  print(DF.rd1000)
}

process.simple.offsets.python.bed<-function(data,prefix="./",window=1000){
  halfw=window/2
  
  pos = round((as.numeric(as.character(data[,2]))+as.numeric(as.character(data[,3])))/2)
  chrs = data[,1]
  chrs = paste("chr",gsub("chr","",chrs),sep="")
  DF.rd1000<-RangedData(ranges=IRanges(start=pos-(halfw-1),end=pos+halfw),space=chrs,ori="",mut="")

  #retrieve reference genome sequences
  DF.wt1000<-as.character(as.matrix(as.data.frame(getSeq(x=Hsapiens,names=as(DF.rd1000,"GRanges")))))
  
  #write to fasta format
  temp<-DNAStringSet(DF.wt1000)
  names(temp)<-paste(DF.rd1000$ori,DF.rd1000$mut,space(DF.rd1000),start(DF.rd1000),end(DF.rd1000),sep = "_")
  writeXStringSet(temp,filepath=paste(prefix,".wt",window,".fasta",sep=""))
  print(DF.rd1000)
}
#read vcf format

if(grepl("vcf$",args[1])){
  data=read.csv(args[1],sep='\t',header=F, comment.char = "#",colClasses = c("character"))
  data[,2]=as.numeric(data[,2])
  data=data[,c(4,5,1,2,2,3)]
  colnames(data)[1:6]<-c("A","B","C","D","E","F")
  data=data[order(data$C,data$D),]
  #for vcf format, we retrieve 1100bp sequences (instead of 1000bp) to allow for small deletions (<100bp)
  process.simple.offsets.python(data,window = 1100,prefix=args[1])
}else if(grepl("bed$",args[1])){
  data=read.csv(args[1],sep='\t',header=F, comment.char = "#",colClasses = c("character"))
  data[,2]=as.numeric(data[,2])
  data[,3]=as.numeric(data[,3])
  colnames(data)[1:3]<-c("A","B","C")
  data=data[order(data$A,data$B),]
  process.simple.offsets.python.bed(data,window = 1000,prefix=args[1])
}
