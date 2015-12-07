#!/usr/bin/env Rscript
library(getopt)
options(digits=3)
Args <- commandArgs()
self<-Args[4]
self<-substring(self,8,nchar(as.character(self)))
spec = c( 
'mask'     , 'm', "0", "character" ," name of brain mask image ",
'input'    , 'i', "1", "character"," the input prefix ", 
'dim'      , 'd', "3", "character"," the input prefix ", 
'output'   , 'o', "output.nii.gz", "character"," the output prefix ")
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
ex<-paste(self," --output myoutput.nii.gz  --input image.nii.gz --mask mask.nii.gz \n \n ")
ex<-format(ex, width=length(ex), justify = c("left"))
cat("\n")
cat(ex)
q(status=1);
}
if ( is.null( opt$dim ) ) dim<-3 else dim<-as.numeric( opt$dim )
#
# take care of optional parameters
for ( myfn in c( opt$mask, opt$input ) )
  {
    if ( !file.exists(myfn) ) 
      {
        print(paste("input file",myfn,"does not exist. Exiting."))
        q(status=1)
      } # else print(paste(" have input file",myfn))
  }
library(ANTsR)
img<-antsImageRead( opt$input, dim )
mask<-antsImageRead( opt$mask, dim )
mymean<-mean( img[ mask == 1] ) 
outimg<-antsImageClone(img)
outimg[ mask == 0 ]<-mymean
outimg = iMath( outimg, "Laplacian", 1, 1 )
antsImageWrite( outimg , opt$output )
