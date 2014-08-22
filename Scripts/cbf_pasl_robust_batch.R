#!/usr/bin/env Rscript
library( "ANTsR" )
library("extremevalues" )
args<-commandArgs(TRUE)

aslImg <- args[1]
m0Img <- args[2]
outname <- args[3]

#init <- cbf_pasl_robust( asl_img, m0_img )
#cbf <- cbf_pasl_quantify( init$delta, init$m0 )

pasl <- antsImageRead( aslImg, 4 )
m0 <- antsImageRead( m0Img, 3 )

pasl.processing <- aslPerfusion(pasl, moreaccurate=FALSE, m0=m0, interpolation="sinc" )
pasl.parameters <- list( sequence="pasl", m0=m0 )
cbf <- quantifyCBF( pasl.processing$perfusion, pasl.processing$mask, pasl.parameters )
meancbf <- cbf$kmeancbf
cbfOverTime <- cbf$timecbf
print(outname)
antsImageWrite( meancbf, outname )

