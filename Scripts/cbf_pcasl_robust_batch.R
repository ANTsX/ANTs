#!/usr/bin/env Rscript
library( "ANTsR" )
library("extremevalues" )
args<-commandArgs(TRUE)

aslImg <- args[1]
outname <- args[2]


pasl <- antsImageRead( aslImg, 4 )

pcasl.processing <- aslPerfusion(pasl, moreaccurate=FALSE)
pcasl.parameters <- list( sequence="pcasl", m0=pcasl.processing$m0 )
cbf <- quantifyCBF( pcasl.processing$perfusion, pcasl.processing$mask, pcasl.parameters )
meancbf <- cbf$kmeancbf
cbfOverTime <- cbf$timecbf
print(outname)
antsImageWrite( meancbf, outname )

