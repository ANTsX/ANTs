#!/usr/bin/env Rscript
library(ANTsR)
library(tools)
if(!usePkg('optparse') | !usePkg('ANTsR')){
  stop("optparse and ANTsR packages required.")
}

optlist <- list(
  make_option(c('-s', '--pCASL'), default='', help=' raw pCASL image'),
  make_option(c('-o', '--outpre'), default='CBF_',
              help='output prefix (defaults to %default)'),
  make_option(c('-t', '--antsCorticalThicknessPrefix'),
              default='', help='prefix of antsCorticalThickness output'),
  make_option(c('-l', '--labelImage'),
              default='', help='label image in template space to warp to ASL'),
  make_option(c('-c', '--paramFile'), default='',
              help='parameter file containing ASL acquisition parameters'),
  make_option(c('-m', '--method'), default='RobustRegression',
              help=paste(' method for perfusion calculation. \n\t\tOne of:',
                '"SimpleSubtraction", "SurroundSubtraction", "SincSubtraction",',
                '"RobustRegression", "BayesianRegression", "LocalBayesianRegression."')),
  make_option(c('-d', '--denoising'), default='CompCorMotion',
                help=paste('denoising method.',
                'Options are: "CompCor", "Motion", "Detrending",',
                '\n\t\t"Cross-Validation", "OutlierRejection".',
                'Multiple options can be specified',
                '(e.g., "CompCorMotion" is legal).  Default is %default.')),
  make_option(c('-b', '--bloodT1'), default=0.67,
              help='blood T1 value (defaults to %default s^-1)'),
  make_option(c('-r', '--robustness'), default=0.95,
              help='robustness parameter (defaults to %default).'),
  make_option(c('-n', '--bootstrapNumber'), default=20,
              help=' number of bootstrap samples (defaults to %default)'),
  make_option(c('-e', '--bootstrapPercent'), default=0.70,
              help='percent to sample per bootstrap run (defaults to %default)'),
  make_option(c('-k', '--keepTmp'), default=F, action='store_true',
              help=paste('keep tmp files, including warps',
                         '(defaults to %default--takes lots of space to save)')),
  make_option(c('-f', '--bootstrapReplace'), default=F, action='store_true',
              help=paste('bootstrap with replacement?  takes arguments',
                         '"false" or "true"; defaults to false.')),
  make_option(c('-v', '--verbose'), default=F, action='store_true',
              help='verbose output.'))

usage <- OptionParser(option_list=optlist, usage='Usage: %prog <s> [apgxetlomdbrnckfv]')
opt <- parse_args(usage)
## debug
#opt <- data.frame(pCASL='/Users/bkandel/CfN//home/bkandel/tmp/pcasl.nii.gz',
#antsCorticalThicknessPrefix=paste(
#      '/Users/bkandel/CfN/data/jag/BD2K01/ASL_pipeline/',
#      'data/AddictionCenter/ABART/',
##      'processed/ABART_Bac_106/Anatomy/ABART_Bac_106_', sep=''),
# labelImage='/Users/bkandel/CfN/data/jag/BD2K01/ASL_pipeline/data/AddictionCenter/template/VS_OFC.nii.gz',
#                  outpre='test')

if(length(opt$verbose > 0)) {
    if(opt$verbose == TRUE) {
    cat('Running antsASLProcessing.R with the following options:\n')
    for(option in names(opt)){
      cat(paste(option, ': ', opt[option], '\n', sep=''))
    }
  }
}

if(length(grep(.Platform$file.sep, opt$outputpre)) > 0) {
  outdir <- dirname(opt$outpre)
  if(!file.exists(outdir)) dir.create(outdir)
}



pcasl <- tryCatch({
    antsImageRead(as.character(opt$pCASL), 4)
  }, error = function(e) {
    stop(paste('pCASL image', as.character(opt$pCASL),
                'does not exist.'))
})

if(file.exists(as.character(opt$paramFile))) {
    config <- read.csv(opt$paramFile)
} else {
    config <- data.frame(tagFirst=T, sequence='pcasl')
}

avg <- getAverageOfTimeSeries(pcasl)
avg <- n3BiasFieldCorrection(avg, 2)
avg <- n3BiasFieldCorrection(avg, 2)
mask <- getMask(avg, mean(avg), Inf, 2)
avg[mask==0] <- 0

moco <- antsMotionCalculation(pcasl, moreaccurate=0)
tag.first <- config$tagFirst
ts <- timeseries2matrix(moco$moco_img, moco$moco_mask)
if (!tag.first) {
  tc <- (rep(c(1, 0), dim(ts)[1])[1:dim(ts)[1]] - 0.5)  # control minus tag
} else {
  tc <- (rep(c(0, 1), dim(ts)[1])[1:dim(ts)[1]] - 0.5)  # tag minus control
}
nuisance <- getASLNoisePredictors(ts, tc)
noise.all <- cbind(moco$moco_params, moco$dvars, nuisance)
noise.combined <- as.matrix(combineNuisancePredictors(ts, tc, noise.all))
censored <- aslCensoring(pcasl, mask, nuis=noise.combined, method='robust')
noise.censored <- noise.combined[censored$which.inliers, ]
perf <- aslAveraging(censored$asl.inlier, mask=moco$moco_mask,
                     nuisance=noise.censored, method='regression')

mvals2 <- apply(ts[tc == 0.5, ], 2, mean)
mvals1 <- apply(ts[tc == -0.5, ], 2, mean)
# mean control should exceed mean tag
if (mean(mvals2) > mean(mvals1)) {
  m0vals<-mvals2
  m1vals<-mvals1
} else {
  m0vals<-mvals1
  m1vals<-mvals2
}

m0 <- antsImageClone(moco$moco_mask)

m0[moco$moco_mask == 0] <- 0
m0[moco$moco_mask == 1] <- m0vals
m0<-n3BiasFieldCorrection(m0,4)
m0<-n3BiasFieldCorrection(m0,2)

cbf <- quantifyCBF(perf, mask=moco$moco_mask,
                   parameters=list(sequence="pcasl", m0=antsImageClone(m0)))
antsImageWrite(cbf$meancbf, paste(opt$outpre, "CBF.nii.gz", sep=""))
antsImageWrite(perf, paste(opt$outpre, "Perfusion.nii.gz", sep=""))

if (nchar(opt$antsCorticalThicknessPrefix) > 0){
  act <- as.character(opt$antsCorticalThicknessPrefix)
  braint1 <- tryCatch({
      antsImageRead(paste(act, "ExtractedBrain0N4.nii.gz", sep=""))
    }, error = function(e) {
      stop(paste('T1 brain image', paste(act, "ExtractedBrain0N4.nii.gz", sep=""),
                  'does not exist.'))
  })
  postnames <- list.files(path=dirname(act),
            glob2rx("*BrainSegmentationPosteriors*.nii.gz"), full.names=TRUE)
  probs <- imageFileNames2ImageList(postnames)
  seg <- tryCatch({
      antsImageRead(paste(act, "BrainSegmentation.nii.gz", sep=""))
    }, error = function(e) {
      stop(paste('Segmentation image', paste(act, "BrainSegmentation.nii.gz", sep=""),
                  'does not exist.'))
  })
  reg.t12asl <- antsRegistration(fixed=avg, moving=braint1,
    typeofTransform="SyNBold", outprefix=as.character(opt$outpre))
  seg.asl <- antsApplyTransforms(avg, seg, reg.t12asl$fwdtransforms, "MultiLabel")
  antsImageWrite(seg.asl, paste(opt$outpre,
                                "SegmentationWarpedToASL.nii.gz", sep=""))
  tx.template2t1 <- c(paste(act, "SubjectToTemplate0GenericAffine.mat", sep=""),
                      paste(act, "SubjectToTemplate1Warp.nii.gz", sep=""))
  tx.template2asl <- c(tx.template2t1, reg.t12asl$invtransforms)
  if (nchar(as.character(opt$labelImage)) > 0) {
    label <- tryCatch( {
      antsImageRead(as.character(opt$labelImage))
    }, error = function(e) {
      stop(paste("Label image", as.character(opt$labelImage),
                 "does not exist."))
    })
    label.asl <- antsApplyTransforms(fixed=avg, moving=label,
      transformlist=tx.template2asl, interpolator="MultiLabel")
    antsImageWrite(label.asl, paste(opt$outpre,
      'LabelWarpedToASL.nii.gz', sep=''))
  }
}
