#!/usr/bin/env Rscript
dotest<-T
options(digits=3)
Args <- commandArgs()
self<-Args[4]
self<-substring(self,8,nchar(as.character(self)))
getPckg <- function(pckg) install.packages(pckg, repos = "http://cran.r-project.org")
pckg = try(require(getopt))
if(!pckg) {
cat("Installing 'getopt' from CRAN\n")
getPckg("getopt")
require("getopt")
}
pckg = try(require(igraph))
if(!pckg) {
  getPckg("igraph")
}
pckg = try(require(glasso))
if(!pckg) {
  getPckg("glasso")
}
library(glasso)
library(igraph)
library(ANTsR)
spec = c( 
'labels'   , 'l', "0", "character" ," name of (aal) label image ", 
'motion'     , 'm', "0", "character" ," name of brain motion csv ",
'mask'     , 'x', "0", "character" ," name of brain mask image ",
'fmri'     , 'f', "0", "character" ," name of BOLD fmri", 
'freq'     , 'q', "0.01x0.1", "character" ," low x high frequency for filtering",
'gdens'    , 'g', 0.25, "numeric","graph density",
'glass'    , 'a', NA, "numeric","graphical lasso parameter",
'help'     , 'h', 0, "logical" ," print the help ", 
'output'   , 'o', "1", "character"," the output prefix ")
# ............................................. #
spec=matrix(spec,ncol=5,byrow=TRUE)
# get the options
opt = getopt(spec)
# ............................................. #
#help was asked for.
if ( !is.null(opt$help) || length(opt) == 1 ) {
#print a friendly message and exit with a non-zero error code
cat("\n")
cat(paste(self,"\n"))
for ( x in 1:nrow(spec) ) {
cat("\n")
longopt<-paste("--",spec[x,1],sep='')
shortopt<-paste("-",spec[x,2],sep='')
hlist<-paste(shortopt,"|",longopt,spec[x,5],"\n \n")
# print(hlist,quote=F)
cat(format(hlist, width=40, justify = c("left")))
}
cat(format("Example: \n", width=40, justify = c("left")))
ex<-paste(self," --output myoutput --mask mask.nii.gz --labels labels.nii.gz --fmri bold.nii.gz  --freq 0.01x0.1   --gdens  0.1 --motion motion.csv \n \n ")
ex<-format(ex, width=length(ex), justify = c("left"))
cat("\n")
cat(ex)
if ( !dotest ) q(status=1);
}
#
for ( myfn in c( opt$mask, opt$fmri, opt$labels , opt$motion ) )
  {
    if ( !file.exists(myfn) ) 
      {
        print(paste("input file",myfn,"does not exist. Exiting."))
        q(status=1)
      } # else print(paste(" have input file",myfn))
  }
freqLo<-0.02
freqHi<-0.1
if ( is.null( opt$gdens ) ) opt$gdens<-0.25
if ( is.null( opt$glass ) ) opt$glass<-NA
if ( ! is.null( opt$freq ) ) {
  freqLo<-as.numeric( strsplit(opt$freq,"x")[[1]][1] )
  freqHi<-as.numeric( strsplit(opt$freq,"x")[[1]][2] )
}
if ( dotest )
  {
    subjid<-"PEDS008_20101120"
    subjid<-"SZ017_20050928"
    subjid<-"NC805_20050906"
    print(paste("start test",subjid,freqHi,freqLo))
    opt$motion<-paste("moco/",subjid,"_MOCOparams.csv",sep='')
    opt$fmri<-paste("bold/",subjid,"_bold.nii.gz",sep='')
    opt$labels<-paste("regions/",subjid,"_BoldRegions.nii.gz",sep='')
    opt$mask<-paste("mask/",subjid,"_Mask.nii.gz",sep='')
    opt$output<-"test/TEST"
    opt$gdens<-0.5
    opt$glass<-0.2
  } else print(paste("start test",opt$output,freqHi,freqLo))
mask<-antsImageRead(opt$mask,3)
aalm<-antsImageRead(opt$labels,3)
bold<-antsImageRead(opt$fmri,4)
motion<-read.csv(opt$motion)
usemotiondirectly<-TRUE 
if ( ! usemotiondirectly ) 
  {
  msvd<-svd( as.matrix(motion[,3:ncol(motion)] ) )
  nsvdcomp <- 5
  motionnuis <- (msvd$u[, 1:nsvdcomp])
  print(paste(" % var of motion ", (sum(msvd$d[1:nsvdcomp])/sum(msvd$d))))
  }
if ( usemotiondirectly ) motionnuis <- as.matrix(motion[,3:ncol(motion)] )
colnames(motionnuis)<-paste("mot",1:ncol(motionnuis),sep='')
negmask<-antsImageClone( mask )
backgroundvoxels <-  negmask == 0
neginds<-which( backgroundvoxels )
negmask[ negmask >= 0 ] <- 0
backgroundvoxels[  ]<-FALSE
backgroundvoxels[ neginds ]<-TRUE
negmask[ backgroundvoxels ]<-1
newnuisnv<-4
bgsvd<-svd(  timeseries2matrix( bold, negmask ) )
mynuis<-bgsvd$u[, 1:newnuisnv]
colnames(mynuis)<-paste("bgdNuis",1:newnuisnv,sep='')
classiccompcor<-compcor(bold,mask=mask,ncompcor = 6)
##################################################
mynuis<-cbind(motionnuis,classiccompcor , mynuis )
print("My nuisance variables are:")
print( colnames(mynuis) )
aalmask<-antsImageClone( aalm )
mylog<-( aalm < 91 & aalm > 0 ) # for aal  
mylog<-( aalm > 0 )
aalmask[ mylog ]<-1
aalmask[!mylog ]<-0
aalm[!mylog]<-0
omat<-timeseries2matrix( bold, aalmask )
mytimes<-dim(omat)[1]
throwaway<-5
motnuisshift<-ashift(motionnuis[throwaway:mytimes,],c(1,0))
motmag<-apply( motionnuis[throwaway:mytimes,], FUN=mean,MARGIN=2)
matmag<-sqrt( sum(motmag[1:9]*motmag[1:9]) )
tranmag<-sqrt( sum(motmag[10:12]*motmag[10:12]) )
motsd<-apply( motionnuis[throwaway:mytimes,]-motnuisshift, FUN=mean,MARGIN=2)
matsd<-sqrt( sum(motsd[1:9]*motsd[1:9]) )
transd<-sqrt( sum(motsd[10:12]*motsd[10:12]) )
omat<-omat[throwaway:mytimes,]
mynuis<-mynuis[throwaway:mytimes,]
mytimes<-dim(omat)[1]
mat<-residuals( lm( omat ~  mynuis ) )
flo<-freqLo
fhi<-freqHi
mytimeshalf<-mytimes/2
# 1st half
mat1<-omat[1:mytimeshalf,]
mynuis1<-mynuis[1:mytimeshalf,]
mat1<-residuals( lm( mat1 ~  mynuis1 ) )
mynetwork1<-filterfMRIforNetworkAnalysis( mat1 , tr=antsGetSpacing(bold)[4], mask=aalmask ,cbfnetwork = "BOLD", labels= aalm , graphdensity = as.numeric(opt$gdens), freqLo = flo, freqHi = fhi , useglasso = opt$glass )
# 2nd half
mat2<-omat[mytimeshalf:mytimes,]
mynuis2<-mynuis[mytimeshalf:mytimes,]
mat2<-residuals( lm( mat2 ~  mynuis2 ) )
mynetwork2<-filterfMRIforNetworkAnalysis( mat2 , tr=antsGetSpacing(bold)[4], mask=aalmask ,cbfnetwork = "BOLD", labels= aalm , graphdensity = as.numeric(opt$gdens), freqLo = flo, freqHi = fhi )
if ( TRUE )
  {
    print( cor.test(mynetwork1$graph$degree,mynetwork2$graph$degree) )
    print( cor.test(mynetwork1$graph$centrality,mynetwork2$graph$centrality) )
    print( cor.test(mynetwork1$graph$closeness,mynetwork2$graph$closeness) )
    print( cor.test(mynetwork1$graph$pagerank,mynetwork2$graph$pagerank) )
    print( cor.test(mynetwork1$graph$betweeness,mynetwork2$graph$betweeness) )
    print( cor.test(mynetwork1$graph$localtransitivity,mynetwork2$graph$localtransitivity) )
  }
cor1<-cor(t(mynetwork1$network))
cor2<-cor(t(mynetwork2$network))
cor1t<-cor1[upper.tri(cor1)]
cor2t<-cor2[upper.tri(cor2)]
# print(cor.test(abs(cor1t),abs(cor2t)))
# print( cor.test(mynetwork1$graph$degree,mynetwork2$graph$degree) )
mynetwork<-filterfMRIforNetworkAnalysis( mat , tr=antsGetSpacing(bold)[4], mask=aalmask ,cbfnetwork = "BOLD", labels= aalm , graphdensity = as.numeric(opt$gdens), freqLo = flo, freqHi = fhi )
if ( FALSE ) {
  par(mfrow=c(1,3))
  pdf(paste(opt$output,'boldrepro.pdf',sep=''),width=5,height=5)
  plot( abs(cor1t),abs(cor2t))
#  plot(mynetwork1$graph$pagerank,mynetwork2$graph$pagerank)
  plot(mynetwork1$graph$degree,mynetwork2$graph$degree)
  plot(mynetwork1$graph$centrality,mynetwork2$graph$centrality)
  dev.off()
}
##################################################
############## now finalize it all ###############
##################################################
reproval<-cor.test(mynetwork1$graph$degree,mynetwork2$graph$degree)$est
names(reproval)<-"boldReproval"
# return network values on full dataset
names(matmag)<-"MatrixMotion"
names(tranmag)<-"TransMotion"
names(matsd)<-"DMatrixMotion"
names(transd)<-"DTransMotion"
mydeg<-mynetwork$graph$degree
names(mydeg)<-paste("Deg",1:length(mynetwork$graph$degree),sep='')
mypr<-mynetwork$graph$pagerank
names(mypr)<-paste("PR",1:length(mynetwork$graph$degree),sep='')
mycent<-mynetwork$graph$centrality
names(mycent)<-paste("Cent",1:length(mynetwork$graph$degree),sep='')
mytrans<-mynetwork$graph$localtransitivity
names(mytrans)<-paste("Trans",1:length(mynetwork$graph$degree),sep='')
myclose<-mynetwork$graph$closeness
names(myclose)<-paste("Close",1:length(mynetwork$graph$degree),sep='')
myc<-c( matmag, tranmag, matsd, transd, reproval , mydeg, mypr, mycent, mytrans, myclose )
outmat<-matrix(myc,nrow=1)
colnames(outmat)<-names(myc)
write.csv(outmat,paste(opt$output,'boldout.csv',sep=''),row.names=F)
write.csv(mynetwork$graph$adjacencyMatrix,paste(opt$output,'adjacencymatrix.csv',sep=''),row.names=F)
write.csv(mynetwork$corrmat,paste(opt$output,'corrmat.csv',sep=''),row.names=F)
write.csv(mynetwork$corrmat,paste(opt$output,'partial_corrmat.csv',sep=''),row.names=F)
if ( dotest ) {
  par(mfrow=c(2,3))
  plot(mynetwork1$graph$degree,mynetwork2$graph$degree)
  plot(mynetwork1$graph$betweeness,mynetwork2$graph$betweeness)
  plot(mynetwork1$graph$closeness,mynetwork2$graph$closeness)
  plot(mynetwork1$graph$pagerank,mynetwork2$graph$pagerank)
  plot(mynetwork2$graph$localtransitivity,mynetwork1$graph$localtransitivity)
  plot(mynetwork1$graph$centrality,mynetwork2$graph$centrality)
}
##################################################
