snp <- read.table("SNP.output.table",header = T)

parents.snp <- snp[(snp$TS.1!=snp$TS.43) &
                  (snp$TS.1=="A/A" | snp$TS.1=="T/T" | snp$TS.1=="C/C" | snp$TS.1=="G/G" |
                    snp$TS.1=="A|A" | snp$TS.1=="T|T" | snp$TS.1=="C|C" | snp$TS.1=="G|G") &
                  (snp$TS.43=="A/A" | snp$TS.43=="T/T" | snp$TS.43=="C/C" | snp$TS.43=="G/G" |
                    snp$TS.43=="A|A" | snp$TS.43=="T|T" | snp$TS.43=="C|C" | snp$TS.43=="G|G") , ]

parents.snp_final <- parents.snp[
  (((parents.snp$TS.1==paste(parents.snp$REF,parents.snp$REF,sep="/")) |
    (parents.snp$TS.1==paste(parents.snp$REF,parents.snp$REF,sep="|")) |
    (parents.snp$TS.1==paste(parents.snp$ALT,parents.snp$ALT,sep="/")) |
    (parents.snp$TS.1==paste(parents.snp$ALT,parents.snp$ALT,sep="|"))) &
    ((parents.snp$TS.43==paste(parents.snp$REF,parents.snp$REF,sep="/")) |
      (parents.snp$TS.43==paste(parents.snp$REF,parents.snp$REF,sep="|")) |
      (parents.snp$TS.43==paste(parents.snp$ALT,parents.snp$ALT,sep="/")) |
      (parents.snp$TS.43==paste(parents.snp$ALT,parents.snp$ALT,sep="|")))),]

chrom <- unique(snp$CHROM)
zhu_198 <- parents.snp_final[,c(1,2,3,4,5)]
zhu_1255 <- parents.snp_final[,c(1,2,3,4,6)]
zhu_1452 <- parents.snp_final[,c(1,2,3,4,7)]
zhu_1746 <- parents.snp_final[,c(1,2,3,4,8)]
TS_43 <- parents.snp_final[,c(1,2,3,4,9)]
TS_1 <- parents.snp_final[,c(1,2,3,4,10)]

heter_homo_draw <- function(special,name){

  k <- 1
  ncolor <- c("slateblue","tomato2")
  mycolors <- c()
  hete_array <- c()
  homo_array <- c()
  a <- 0

  for (i in chrom){
    chromm <- special[special$CHROM==i,]
    j <- 0
    WindowLength <- 1000000

    while(j <= chromm$POS[length(chromm$POS)]){
      eWindow <- chromm[chromm$POS >= j & chromm$POS < j + WindowLength, ]

      #hete <- sum(eWindow == paste(eWindow$REF, eWindow$ALT,sep="/") |
      #              eWindow == paste(eWindow$REF, eWindow$ALT,sep="|") |
      #              eWindow == paste(eWindow$ALT, eWindow$REF,sep="/") |
      #              eWindow == paste(eWindow$ALT, eWindow$REF,sep="|")) / WindowLength

      hete <- sum((eWindow[,5] != paste(eWindow$REF, eWindow$REF,sep="/")) &
                    (eWindow[,5] != paste(eWindow$REF, eWindow$REF,sep="|")) &
                    (eWindow[,5] != paste(eWindow$ALT, eWindow$ALT,sep="/")) &
                    (eWindow[,5] != paste(eWindow$ALT, eWindow$ALT,sep="|"))) / length(eWindow[,5])

      homo <- sum((eWindow[,5] == paste(eWindow$REF, eWindow$REF,sep="/")) |
                    (eWindow[,5] == paste(eWindow$REF, eWindow$REF,sep="|")) |
                    (eWindow[,5] == paste(eWindow$ALT, eWindow$ALT,sep="/")) |
                    (eWindow[,5] == paste(eWindow$ALT, eWindow$ALT,sep="|"))) / length(eWindow[,5])

      hete_array <- c(hete_array, hete)
      homo_array <- c(homo_array, homo)
      mycolors <- c(mycolors,ncolor[k])
      j <- j + WindowLength
    }
    print(hete_array)

    a <- c(a, length(hete_array))
    if (k == 1){
      k <- 2
    }
    else{
      k <- 1
    }
  }

  p <- c()
  for (i in c(2:14)){
    p <- c(p,a[i-1]+(a[i]-a[i-1])/2)
  }
  maxHete_lim <- max(hete_array)
  maxHomo_lim <- max(homo_array)
  if (name != "TS_1" & name != "TS_43"){
    i <- 0
    while(maxHete_lim <= 1){
      i <- i + 1
      maxHete_lim <- max(hete_array) * 10^i
    }
  }
  j <- 0
  while(maxHomo_lim <= 1){
    j <- j + 1
    maxHomo_lim <- max(homo_array) * 10^j
  }
  pdf(paste(name,".pdf"),width = 9,height = 7)
  par(mfrow=c(2,1))
  plot(hete_array,cex=0.5,ylab = ("Heterozygosity"),cex.axis=0.7,
       cex.lab=0.7,xlab = "Chromosome(Mb)", col=mycolors,xaxt="n",
       ylim=c(0,1),main=paste("Heterozygosity of", name,sep = ' '))
  axis(side=1,at=p,cex.axis=0.7,labels=c("chr00","chr01","chr02","chr03","chr04","chr05",
                                         "chr06","chr07","chr08","chr09","chr10","chr11","chr12"))
  #axis(side=2,at=seq(0,ceiling(maxHete_lim)/10**i,ceiling(maxHete_lim)/10**i/4),
  #     cex.axis=0.5,labels=seq(0,ceiling(maxHete_lim),ceiling(maxHete_lim)/4))

  plot(homo_array,cex=0.5,ylab = ("Similarity"),cex.axis=0.7,
       cex.lab=0.7,xlab = "Chromosome(Mb)", col=mycolors,xaxt="n",
       ylim=c(0,1),main=paste("Similarity of", name,sep = ' '))
  axis(side=1,at=p,cex.axis=0.7,labels=c("chr00","chr01","chr02","chr03","chr04","chr05",
                                         "chr06","chr07","chr08","chr09","chr10","chr11","chr12"))
  #axis(side=2,at=seq(0,ceiling(maxHomo_lim)/10**j,ceiling(maxHomo_lim)/10**j/4),
  #     cex.axis=0.5,labels=seq(0,ceiling(maxHomo_lim),ceiling(maxHomo_lim)/4))
  dev.off()
}

heter_homo_draw(zhu_198,"zhu_198")
heter_homo_draw(zhu_1746,"zhu_1746")
heter_homo_draw(zhu_1255,"zhu_1255")
heter_homo_draw(zhu_1452,"zhu_1452")
heter_homo_draw(TS_43,"TS_43")
heter_homo_draw(TS_1,"TS_1")
