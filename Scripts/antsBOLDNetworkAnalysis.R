#!/usr/bin/env Rscript
myscale <- function(x,doscale=F) {
    if ( doscale ) return( scale(x) )
    return( x )
}
dotest<-F
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
pckg = try(require(psych))
if(!pckg) {
  getPckg("psych")
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
'thresh'   , 't', "1x100000000", "character" ," lower upper thresh for label image ",
'motion'     , 'm', "0", "character" ," name of brain motion csv ",
'mask'     , 'x', "0", "character" ," name of brain mask image ",
'fmri'     , 'f', "0", "character" ," name of BOLD fmri",
'freq'     , 'q', "0.01x0.1", "character" ," low x high frequency for filtering",
'gdens'    , 'g', 0.25, "numeric","graph density",
'winsortrim' , 'w', 0.01, "numeric","winsorizing value e.g. 0.05 = 5%",
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
threshLo<-1
threshHi<-10000000
throwaway<-10
if ( is.null( opt$winsortrim ) ) opt$winsortrim<-0
if ( is.null( opt$gdens ) ) opt$gdens<-0.25
if ( is.null( opt$glass ) ) opt$glass<-NA
if ( ! is.null( opt$freq ) ) {
  freqLo<-as.numeric( strsplit(opt$freq,"x")[[1]][1] )
  freqHi<-as.numeric( strsplit(opt$freq,"x")[[1]][2] )
}
if ( ! is.null( opt$thresh ) ) {
  threshLo<-as.numeric( strsplit(opt$thresh,"x")[[1]][1] )
  threshHi<-as.numeric( strsplit(opt$thresh,"x")[[1]][2] )
}
if ( dotest )
  {
    subjid<-"SZ017_20050928"
    subjid<-"NC805_20050906"
#    subjid<-"PEDS045_20110208"
#    subjid<-"PEDS008_20101120"
#    subjid<-"PEDS104_20120817"
    subjid<-"PEDS007_20110903"
    print(paste("start test",subjid,freqHi,freqLo))
    opt$motion<-paste("moco/",subjid,"_MOCOparams.csv",sep='')
    opt$fmri<-paste("bold/",subjid,"_bold.nii.gz",sep='')
    opt$labels<-paste("regions/",subjid,"_AAL.nii.gz",sep='')
    opt$mask<-paste("mask/",subjid,"_Mask.nii.gz",sep='')
    opt$output<-"test/TEST"
    opt$gdens<-0.25
    opt$glass<-0.01
    threshLo<-1
    threshHi<-90
    freqLo<-0.01
    freqHi<-0.1
  } else print(paste("start subject",opt$output,freqHi,freqLo))
mask<-antsImageRead(opt$mask,3)
aalm<-antsImageRead(opt$labels,3)
bold<-antsImageRead(opt$fmri,4)
mytimes<-dim(bold)[4]
motion<-read.csv(opt$motion)
aalmask<-antsImageClone( aalm )
mylog<-( aalm >= threshLo & aalm <= threshHi )
aalmask[ mylog ]<-1
aalmask[!mylog ]<-0
aalm[!mylog]<-0
print(paste("You are using:",length( unique( aalm[aalmask>0] ) ) ,"unique labels."))
omat<-myscale( timeseries2matrix( bold, aalmask ) )
templateFD<-rep(0,nrow(motion))
DVARS<-rep(0,nrow(motion))
for ( i in 2:nrow(motion) ) {
  mparams1<-c( motion[i,3:14] )
  tmat1<-matrix( as.numeric(mparams1[1:9]), ncol = 3, nrow = 3)
  mparams2<-c( motion[i-1,3:14] )
  tmat2<-matrix( as.numeric(mparams2[1:9]), ncol = 3, nrow = 3)
  pt<-t( matrix(  rep(10,3), nrow=1) )
  newpt1<-data.matrix(tmat1) %*%  data.matrix( pt )+as.numeric(mparams1[10:12])
  newpt2<-data.matrix(tmat2) %*%  data.matrix( pt )+as.numeric(mparams1[10:12])
  templateFD[i]<-sum(abs(newpt2-newpt1))
  DVARS[i]<-sqrt( mean( ( omat[i,] - omat[i-1,] )^2 ) )
}
# question - should this be a constant of 0.2 as recommended in Yan Craddock He Milham?
keepinds<-which( templateFD < ( mean(templateFD) + 2*sd(templateFD)) & ( (1:mytimes) > throwaway ) )
keepinds<-c(throwaway,keepinds)
throwinds<-which( templateFD > ( mean(templateFD) + 2*sd(templateFD)) & ( (1:mytimes) > throwaway ) )
if ( dotest ) { plot( templateFD , type='l' ); print(paste(sum( templateFD > 0.2 ),length(throwinds))) }
doimpute<-TRUE
if ( length( throwinds )  > 0 & doimpute )
for ( i in throwinds ) {
  previ <- max( keepinds[ keepinds < i ] )
  nexti <- min( keepinds[ keepinds > i ] )
  wt1 <-  1-abs( i - previ )/(nexti-previ)
  wt2 <-  1-abs( i - nexti )/(nexti-previ)
  omat[i,] <- wt1 * omat[previ,] + omat[nexti,] * wt2
  DVARS[i] <- wt1 * DVARS[previ] + DVARS[nexti] * wt2
  templateFD[i] <- wt1 * templateFD[previ] + templateFD[nexti] * wt2
  motion[i,] <- wt1 * motion[previ,] + motion[nexti,] * wt2
  }
keepinds<-throwaway:mytimes
usemotiondirectly<-TRUE
if ( ! usemotiondirectly )
  {
  msvd<-svd( as.matrix(motion[keepinds,3:ncol(motion)] ) )
  mysum<-cumsum(msvd$d)/sum(msvd$d)
  nsvdcomp<-which( mysum > 0.999 )[1]
  motionnuis <- (msvd$u[, 1:nsvdcomp])
  print(paste(" % var of motion ", mysum[nsvdcomp],'with',nsvdcomp ) )
  }
if ( usemotiondirectly ) motionnuis <- as.matrix(motion[keepinds,3:ncol(motion)] )
colnames(motionnuis)<-paste("mot",1:ncol(motionnuis),sep='')
bkgd<-4
if ( bkgd > 0 ) {
  negmask<-antsImageClone( mask )
  backgroundvoxels <-  negmask == 0
  neginds<-which( backgroundvoxels )
  negmask[ negmask >= 0 ] <- 0
  backgroundvoxels[  ]<-FALSE
  backgroundvoxels[ neginds ]<-TRUE
  negmask[ backgroundvoxels ]<-1
  ImageMath(3,negmask,"ME",negmask,1)
  tempmat<-myscale( timeseries2matrix( bold, negmask )[keepinds,] )
#  tempmat<-myscale( timeseries2matrix( bold, mask )[keepinds,] )
#  tempmat<-tempmat-ashift(tempmat,c(1,0))
  bgsvd<-svd( tempmat )
  mysum<-cumsum(bgsvd$d)/sum(bgsvd$d)
  newnuisv<-min( c( bkgd, which( mysum > 0.8 )[1] ) )
  print(paste(newnuisv," % var of bgd ",mysum[newnuisv] ) )
  bgdnuis<-bgsvd$u[, 1:newnuisv]
  colnames(bgdnuis)<-paste("bgdNuis",1:newnuisv,sep='')
}
print(paste("winsorizing with trim",opt$winsortrim))
if ( opt$winsortrim > 0 ) omat<-winsor(omat,trim=opt$winsortrim)
omat<-omat[keepinds,]
##################################################
classiccompcor<-compcor(omat,mask=mask,ncompcor = 4 )
omotionnuis<-as.matrix(motion[keepinds,3:ncol(motion)] )
motnuisshift<-ashift(omotionnuis,c(1,0))
motmag<-apply( omotionnuis, FUN=mean,MARGIN=2)
matmag<-sqrt( sum(motmag[1:9]*motmag[1:9]) )
tranmag<-sqrt( sum(motmag[10:12]*motmag[10:12]) )
motsd<-apply( omotionnuis-motnuisshift, FUN=mean,MARGIN=2)
matsd<-sqrt( sum(motsd[1:9]*motsd[1:9]) )
transd<-sqrt( sum(motsd[10:12]*motsd[10:12]) )
dmatrix<-(omotionnuis-motnuisshift)[,1:9]
dtran<-(omotionnuis-motnuisshift)[,10:12]
dmatrixm<-apply( dmatrix * dmatrix , FUN=sum, MARGIN=1 )
dtranm<-apply( dtran * dtran , FUN=sum, MARGIN=1 )
if ( bkgd )
mynuis<-cbind(scale(dmatrixm)[,1],scale(dtranm)[,1],classiccompcor, bgdnuis, templateFD[keepinds], DVARS[keepinds] ) else mynuis<-cbind(scale(dmatrixm)[,1],scale(dtranm)[,1],classiccompcor, templateFD[keepinds], DVARS[keepinds] )
colnames(mynuis)[1:2]<-c("dmatrix","dtran")
colnames(mynuis)[(length(colnames(mynuis))-1):length(colnames(mynuis))]<-c("FD","DVARS")
print("My nuisance variables are:")
print( colnames(mynuis) )
mytimes<-dim(omat)[1]
mat<-myscale( residuals( lm( omat ~  mynuis ) ) , doscale = TRUE )
flo<-freqLo
fhi<-freqHi
mytimeshalf<-mytimes/2
# 1st half
mat1<-omat[1:mytimeshalf,]
mynuis1<-mynuis[1:mytimeshalf,]
mat1<-residuals( lm( mat1 ~  mynuis1 ) )
locmotnuis<-cbind( mynuis, motionnuis-ashift(motionnuis,c(1,0)) )
mynetwork1<-filterfMRIforNetworkAnalysis( mat1 , tr=antsGetSpacing(bold)[4], mask=aalmask ,cbfnetwork = "BOLD", labels= aalm , graphdensity = as.numeric(opt$gdens), freqLo = flo, freqHi = fhi , useglasso = opt$glass, usesvd=FALSE  ) # , nuisancein = locmotnuis[1:mytimeshalf,])
# 2nd half
mat2<-omat[mytimeshalf:mytimes,]
mynuis2<-mynuis[mytimeshalf:mytimes,]
mat2<-residuals( lm( mat2 ~  mynuis2 ) )
mynetwork2<-filterfMRIforNetworkAnalysis( mat2 , tr=antsGetSpacing(bold)[4], mask=aalmask ,cbfnetwork = "BOLD", labels= aalm , graphdensity = as.numeric(opt$gdens), freqLo = flo, freqHi = fhi, usesvd=FALSE   ) # , nuisancein = locmotnuis[mytimeshalf:mytimes,] )
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
mynetwork<-filterfMRIforNetworkAnalysis( mat , tr=antsGetSpacing(bold)[4], mask=aalmask ,cbfnetwork = "BOLD", labels= aalm , graphdensity = as.numeric(opt$gdens), freqLo = flo, freqHi = fhi , usesvd=FALSE ) # , nuisancein = locmotnuis )
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
reproval2<-cor.test(mynetwork1$graph$localtransitivity,mynetwork2$graph$localtransitivity)$est
names(reproval)<-"boldReproval"
names(reproval2)<-"boldReproval2"
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
mybtwn<-mynetwork$graph$betweeness
names(mybtwn)<-paste("Btwn",1:length(mynetwork$graph$degree),sep='')
mytimes<-length(keepinds)
names(mytimes)<-"NTimePoints"
meanFD<-mean(templateFD)
names(meanFD)<-"meanFD"
meanDVARS<-mean(DVARS)
names(meanDVARS)<-"meanDVARS"
myc<-c( matmag, tranmag, matsd, transd, reproval , mydeg, mypr, mycent, mytrans, myclose , mybtwn, reproval2, mytimes, meanFD, meanDVARS )
outmat<-matrix(myc,nrow=1)
colnames(outmat)<-names(myc)
write.csv(outmat,paste(opt$output,'boldout.csv',sep=''),row.names=F)
write.csv(mynetwork$graph$adjacencyMatrix,paste(opt$output,'adjacencymatrix.csv',sep=''),row.names=F)
write.csv(mynetwork$corrmat,paste(opt$output,'pearson_corrmat.csv',sep=''),row.names=F)
write.csv(mynetwork$partialcorrmat,paste(opt$output,'partial_corrmat.csv',sep=''),row.names=F)
write.csv(mynetwork$glassocormat,paste(opt$output,'glasso_corrmat.csv',sep=''),row.names=F)
write.csv(mynetwork$rcormat,paste(opt$output,'r_corrmat.csv',sep=''),row.names=F)
if ( TRUE ) {
  require(gridExtra)
  require(ggplot2)
  dat<-data.frame(deg1=mynetwork1$graph$degree,deg2=mynetwork2$graph$degree)
  plot1<-ggplot(dat, aes(x=deg1, y=deg2)) + geom_point(shape=1) + geom_smooth(method=lm)
  dat<-data.frame(btwn1=mynetwork1$graph$betweeness,btwn2=mynetwork2$graph$betweeness)
  plot2<-ggplot(dat, aes(x=btwn1, y=btwn2)) + geom_point(shape=1) + geom_smooth(method=lm)
  dat<-data.frame(clos1=mynetwork1$graph$closeness,clos2=mynetwork2$graph$closeness)
  plot3<-ggplot(dat, aes(x=clos1, y=clos2)) + geom_point(shape=1) + geom_smooth(method=lm)
  dat<-data.frame(pr1=mynetwork1$graph$pagerank,pr2=mynetwork2$graph$pagerank)
  plot4<-ggplot(dat, aes(x=pr1, y=pr2)) + geom_point(shape=1) + geom_smooth(method=lm)
  dat<-data.frame(LT1=mynetwork1$graph$localtransitivity,LT2=mynetwork2$graph$localtransitivity)
  plot5<-ggplot(dat, aes(x=LT1, y=LT2)) + geom_point(shape=1) + geom_smooth(method=lm)
  dat<-data.frame(cent1=mynetwork1$graph$centrality,cent2=mynetwork2$graph$centrality)
  plot6<-ggplot(dat, aes(x=cent1, y=cent2)) + geom_point(shape=1) + geom_smooth(method=lm)
  pdf(paste(opt$output,"reproplot.pdf",sep=''),width=8,height=4)
  grid.arrange(plot1, plot2,plot3, plot4,plot5, plot6, ncol=3)
  dev.off()
}
##################################################
