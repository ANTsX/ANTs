#!/usr/bin/env Rscript
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
library(igraph)
#
# helper functions
#
avgimg <- function(  mylist , mask )
{
avg<-antsImageClone( mylist[[1]] )
avg[ mask == 1 ]<-0
for ( i in 1:length(mylist) ) 
  {
  avg[ mask == 1 ] <- avg[ mask == 1 ] + mylist[[i]][ mask == 1 ] * 1/length(mylist)
  }
 return( avg )
}
#
sdimg <- function(  mylist , mask )
{
avg<-avgimg( mylist , mask )
sdi<-antsImageClone( avg )
sdi[ mask == 1 ]<-0
for ( i in 1:length(mylist) ) 
  {
  sdi[ mask == 1 ] <- sdi[ mask == 1 ] + abs( mylist[[i]][ mask == 1 ] - avg[ mask == 1 ] )  * 1/length(mylist)
  }
 return( sdi )
}
#
interleave <- function(v1,v2)
{
ord1 <- 2*(1:length(v1))-1
ord2 <- 2*(1:length(v2))
c(v1,v2)[order(c(ord1,ord2))]
}
#
#
spec = c( 
'labels'   , 'l', "0", "character" ," name of (aal) label image ", 
'mask'     , 'x', "0", "character" ," name of brain mask image ",
'seg'      , 's', "0", "character" ," name of 3-tissue segmentation image 2 == gm ", 
'fmri'     , 'f', "0", "character" ," name of fmri (ASL or BOLD)", 
'modality' , 'w', "BOLD", "character" ," which modality ASLCBF, ASLBOLD or BOLD ",
'freq'     , 'q', "0.01x0.1", "character" ," low x high frequency for filtering",
'gdens'    , 'g', 0.25, "numeric","graph density",
'tr'       , 't', "2x4", "character","TR for BOLD and ASL e.g. 2.2x4",
'help'     , 'h', 0, "logical" ," print the help ", 
'output'   , 'o', "1", "character"," the output prefix ", 
'bloodt1'  , 'b', 2, "numeric", "blood relaxation (inv of t1, defaults to 0.67 s^-1", 
'robust'   , 'r', 2, "numeric", "robustness parameter", 
'nboot'    , 'n', 2, "numeric", "number of bootstrap runs", 
'pctboot'  , 'p', 2, "numeric", "percent to sample per bootstrap run", 
'replace ' , 'e', 2, "logical", "resample with replacement during bootstrap?")
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
ex<-paste(self," -o myoutput -x mask.nii.gz --labels labels.nii.gz --fmri bold.nii.gz --modality ASLCBF  --freq 0.03x0.08 \n \n ")
ex<-format(ex, width=length(ex), justify = c("left"))
cat("\n")
cat(ex)
q(status=1);
}
#
# take care of optional parameters
if(is.null(opt$bloodt1)) opt$bloodt1 <- 0.67
if(is.null(opt$robust)) opt$robust <- 0.95
if(is.null(opt$nboot)) opt$nboot <- 20
if(is.null(opt$pctboot)) opt$pctboot <- 0.70
if(opt$pctboot > 1.0) {
  cat('pctboot was greater than 1; setting to 70%.\n')
  opt$pctboot <- 0.70 
}
if(is.null(opt$replace)) opt$replace <- FALSE
for ( myfn in c( opt$mask, opt$fmri, opt$labels ) )
  {
    if ( !file.exists(myfn) ) 
      {
        print(paste("input file",myfn,"does not exist. Exiting."))
        q(status=1)
      } # else print(paste(" have input file",myfn))
  }
freqLo<-0.01
freqHi<-0.1
trASL<-4
trBOLD<-2.2
if ( is.null( opt$gdens ) ) opt$gdens<-0.25
if ( ! is.null( opt$freq ) ) {
  freqLo<-as.numeric( strsplit(opt$freq,"x")[[1]][1] )
  freqHi<-as.numeric( strsplit(opt$freq,"x")[[1]][2] )
}
if ( ! is.null( opt$tr ) ) {
  trASL<-as.numeric( strsplit(opt$tr,"x")[[1]][2] )
  trBOLD<-as.numeric( strsplit(opt$tr,"x")[[1]][1] )
}
print( paste( "fLo",freqLo,"fHi", freqHi ,"modality",opt$modality,"graphdensity", opt$gdens ) )
suppressMessages( library(ANTsR) )
seg<-NA
if ( ! is.null( opt$seg ) ) seg<-antsImageRead( opt$seg , 3 )
fmri<-antsImageRead( opt$fmri, 4 )
aal2fmri<-antsImageRead( opt$labels, 3 )
mask<-antsImageRead( opt$mask, 3 )
if ( as.character(opt$modality) == "ASLCBF" | as.character(opt$modality) == "ASLBOLD" )
  {
    mat<-timeseries2matrix( fmri, mask )
    cbflist<-list( )
    moco_results <- motion_correction(fmri)
    regweights <- aslPerfusion(fmri, mask=mask, moreaccurate=T, dorobust=opt$robust, moco_results=moco_results)$regweights
    for ( i in 1:opt$nboot ) {
      timeinds<-sample( 2:nrow(mat) , round( nrow(mat) )*(opt$pctboot/2) , replace=opt$replace ) 
      timeinds<-( timeinds %% 2 )+timeinds
      timeinds<-interleave( timeinds-1, timeinds )
      aslarr<-as.array( fmri ) 
      aslarr2<-aslarr[,,,timeinds]
      aslsub<-as.antsImage( aslarr2 )
      antsSetSpacing( aslsub , antsGetSpacing( fmri ) )
      
      mocoarr <- as.array(moco_results$moco_img)
      mocoarr2<-mocoarr[,,,timeinds]
      mocosub<-as.antsImage( mocoarr2 )
      antsSetSpacing( mocosub , antsGetSpacing( fmri ) )

      mocoparams <- as.data.frame(moco_results$moco_params)
      mocoparams.sub <- mocoparams[timeinds, ]

      moco_results.sub <- list(moco_img=mocosub, moco_params=mocoparams.sub, moco_avg_img=moco_results$moco_avg_img)
      regweights.sub <- regweights[timeinds] 
      proc <- aslPerfusion( aslsub, mask=mask, moreaccurate=TRUE, dorobust=opt$robust, moco_results=moco_results.sub, regweights=regweights.sub)
      param <- list( sequence="pcasl", m0=proc$m0 )
      cbf <- quantifyCBF( proc$perfusion, mask, param )
      cbflist<-lappend( cbflist, cbf$kmeancbf )
    }
    write.csv(data.frame(Regweights=regweights), paste(opt$output, 'ExcludedTimePoints.csv', sep=''))
    motion <- as.data.frame(moco_results$moco_params)
    templateFD<-rep(0,nrow(motion))
    DVARS<-rep(0,nrow(motion))
    omat <- mat 
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
    omotionnuis<-as.matrix(motion[, 3:ncol(motion)] )
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
    names(matmag)<-"MatrixMotion"
    names(tranmag)<-"TransMotion"
    names(matsd)<-"DMatrixMotion"
    names(transd)<-"DTransMotion"
    write.csv(cbind(motion, data.frame(templateFD=templateFD, DVARS=DVARS)), paste(opt$output, 'MotionParams.csv', sep=''))
    write.csv(data.frame(MatrixMotion=matmag, Transmotion=tranmag, DMatrixMotion=matsd, DTransMotion=transd), 
      paste(opt$output, 'MotionSummary.csv', sep=''), row.names=F)
    cbfout<-antsImageClone( mask )
    avgcbf<-avgimg( cbflist , mask )
    sdi<-sdimg( cbflist , mask )
    thresh <- 75
    cbfout[ sdi > thresh ] <- 0
    cbfout[ sdi <= thresh ] <- avgcbf[ sdi <= thresh ]
    fn<-paste( opt$output,"_kcbf.nii.gz",sep='')
    antsImageWrite( cbfout , fn )
    pcasl.processing <- aslPerfusion( fmri, mask=mask, moreaccurate=TRUE , dorobust = opt$robust, moco_results=moco_results)
    pcasl.parameters <- list( sequence="pcasl", m0=pcasl.processing$m0, T1b=opt$bloodt1)
    cbf <- quantifyCBF( pcasl.processing$perfusion, mask, pcasl.parameters )
    filterpcasl<-getfMRInuisanceVariables( fmri, mask = mask , moreaccurate=TRUE )
    xideal<-pcasl.processing$xideal
    tsResid<-residuals( lm( filterpcasl$matrixTimeSeries ~ filterpcasl$nuisancevariables + xideal ))
    mynetwork<-filterfMRIforNetworkAnalysis( tsResid , tr=trASL, mask=mask ,cbfnetwork = opt$modality, labels = aal2fmri , graphdensity = opt$gdens, freqLo = freqLo, freqHi = freqHi , seg = seg )
  }
if ( as.character(opt$modality) == "BOLD" )
  {
    dd<-getfMRInuisanceVariables( fmri, mask = mask , moreaccurate=TRUE )
    ###### do a simple thing instead
    mask<-dd$mask
    negmask<-antsImageClone( mask )
    backgroundvoxels <-  negmask == 0
    neginds<-which( backgroundvoxels )
    neginds<-sample(neginds,length(neginds)/25)
    negmask[ negmask >= 0 ] <- 0
    backgroundvoxels[  ]<-FALSE
    backgroundvoxels[ neginds ]<-TRUE
    negmask[ backgroundvoxels ]<-1
    mynuis<-svd(  timeseries2matrix( fmri, negmask ) )$u[, 1:8]
    colnames(mynuis)<-paste("bgdNuis",1:8,sep='')
    dd$nuisancevariables<-cbind(dd$nuisancevariables[,1:4],mynuis)
    
    tsResid<-residuals( lm( dd$matrixTimeSeries ~ dd$nuisancevariables  ))  # see http://www.ncbi.nlm.nih.gov/pubmed/21889994
    mynetwork<-filterfMRIforNetworkAnalysis( tsResid , tr=trBOLD, mask=dd$mask ,cbfnetwork = "BOLD", labels = aal2fmri , graphdensity = opt$gdens, freqLo = freqLo, freqHi = freqHi , seg = seg )
  }
print( names( mynetwork ) )
fn<-paste( opt$output,opt$modality,"_corrmat.mha",sep='')
print( paste( "write correlation matrix",fn ) )
antsImageWrite( as.antsImage( mynetwork$corrmat ) , fn )
g<-mynetwork$graph
print( names( g ) )
fn<-paste( opt$output,opt$modality,"_graph_metrics.csv",sep='')
print( paste( "write graph metrics",fn ) )
graphmetrics<-data.frame( degree =      g$degree   ,           closeness = g$closeness ,
                          pagerank =    g$pagerank ,         betweenness = g$betweeness ,
                          localtransitivity = g$localtransitivity ,
                          modularity        = g$walktrapcomm$modularity )
write.csv( graphmetrics , fn, row.names = F , quote = F )
 
