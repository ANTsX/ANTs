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
ex<-paste(self," -o myoutput -x mask.nii.gz --labels labels.nii.gz --fmri bold.nii.gz --modality ASLCBF  --freq 0.03x0.08 \n \n ")
ex<-format(ex, width=length(ex), justify = c("left"))
cat("\n")
cat(ex)
q(status=1);
}
#
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
    pcasl.processing <- aslPerfusion( fmri, mask=mask, moreaccurate=TRUE , dorobust = 0.85 )
    pcasl.parameters <- list( sequence="pcasl", m0=pcasl.processing$m0 )
    cbf <- quantifyCBF( pcasl.processing$perfusion, mask, pcasl.parameters )
    fn<-paste( opt$output,"_kcbf.nii.gz",sep='')
    antsImageWrite( cbf$kmeancbf , fn )
    filterpcasl<-getfMRInuisanceVariables( fmri, mask = mask , moreaccurate=TRUE )
    xideal<-pcasl.processing$xideal
    tsResid<-residuals( lm( filterpcasl$matrixTimeSeries ~ filterpcasl$nuisancevariables + xideal ))
    mynetwork<-filterfMRIforNetworkAnalysis( tsResid , tr=trASL, mask=mask ,cbfnetwork = opt$modality, labels = aal2fmri , graphdensity = opt$gdens, freqLo = freqLo, freqHi = freqHi , seg = seg )
  }
if ( as.character(opt$modality) == "BOLD" )
  {
    dd<-getfMRInuisanceVariables( fmri, mask = mask , moreaccurate=TRUE )
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
 
