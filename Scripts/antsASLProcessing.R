#!/usr/bin/env Rscript
library(ANTsR)
library(tools)
if(!usePkg('optparse') | !usePkg('ANTsR')){
  stop("optparse and ANTsR packages required.")
}

optlist <- list(
  make_option(c('-s', '--pCASL'), default='', help=' raw pCASL image'),
  make_option(c('-o', '--outputpre'), default='CBF_',
              help='output prefix (defaults to %default)'),
  make_option(c('-a', '--antsCorticalThicknessPrefix'),
              default='', help='prefix of antsCorticalThickness output'),
  make_option(c('-l', '--labelSet'),
              default='', help='label set in template space to warp to ASL'),
  make_option(c('-t', '--template'),
              default='', help='Template to warp output to'),
  make_option(c('-c', '--paramFile'), default='',
              help='parameter file containing ASL acquisition parameters'),
  make_option(c('-x', '--smoothingFWHM'), default=0,
              help='Full width half max for smoothing'),
  make_option(c('-m', '--method'), default='regression',
              help=paste(' method for perfusion calculation. \n\t\tOne of:',
                '"regression", "subtraction", "bayesian",',
                '"RobustRegression", "BayesianRegression", "LocalBayesianRegression."')),
  make_option(c('-d', '--denoising'), default='CompCorMotion',
                help=paste('denoising method.',
                'Options are: "CompCor", "Motion", "Detrending",',
                '\n\t\t"Cross-Validation", "OutlierRejection".',
                'Multiple options can be specified',
                '(e.g., "CompCorMotion" is legal).  Default is %default.')),
  make_option(c('-g', '--debug'), default=0,
               help=paste('Save debugging information, including motion',
                  'correction and nuisance variables')),
  make_option(c('-b', '--bloodT1'), default=0.67,
              help='blood T1 value (defaults to %default s^-1)'),
  make_option(c('-r', '--robustness'), default=0.95,
              help='robustness parameter (defaults to %default).'),
  make_option(c('-n', '--bootstrapNumber'), default=20,
              help=' number of bootstrap samples (defaults to %default)'),
  make_option(c('-e', '--bootstrapPercent'), default=0.70,
              help='percent to sample per bootstrap run (defaults to %default)'),
  make_option(c('-k', '--keepTmp'), default=0,
              help=paste('keep tmp files, including warps',
                         '(defaults to %default--takes lots of space to save)')),
  make_option(c('-f', '--bootstrapReplace'), default=0,
              help=paste('bootstrap with replacement?  takes arguments',
                         '0 or 1; defaults to 0.')),
  make_option(c('-v', '--verbose'), default=0,
              help='verbose output.'))

usage <- OptionParser(option_list=optlist, usage='Usage: %prog <s> [otlcxmdgbrnekfv]')
opt <- parse_args(usage)
## debug
#opt <- data.frame(
# pCASL=paste('/data/jag/BD2K01/ASL_pipeline/data/AddictionCenter/ABART/imgs/',
#             '../processed/ABART_Bac_106/ASL/ABART_Bac_106_pCASL.nii.gz', sep=''),
# outputpre=paste('/data/jag/BD2K01/ASL_pipeline/data/AddictionCenter/ABART/imgs',
#    '/../processed/ABART_Bac_106/ASL/ABART_Bac_106_', sep=''),
# antsCorticalThicknessPrefix=paste('/data/jag/BD2K01/ASL_pipeline/',
#      'data/AddictionCenter/ABART/imgs/../processed/ABART_Bac_106',
#      '/ASL/../Anatomy/ABART_Bac_106_', sep=''),
# labelSet=paste('/data/jag/BD2K01/ASL_pipeline/templates/',
#                'HarvardOxford/ABART_rois.nii.gz', sep=''),
# template=paste('/data/jag/BD2K01/ASL_pipeline/templates/',
#                'HarvardOxford/MNI152_T1_2mm.nii.gz', sep=''))
#                  pCASL='data/101_pcasl.nii.gz',
#                  out='test')


if(!file.exists(as.character(opt$pCASL))) {
  stop(paste('pCASL image', opt$pCASL,
    'does not exist.'))
}

if(opt$verbose) {
  cat('Running antsASLProcessing.R with the following options:\n')
  for(option in names(opt)){
    cat(paste(option, ': ', opt[option], '\n', sep=''))
  }
}

if(length(grep(.Platform$file.sep, opt$outputpre)) > 0) {
  outdir <- dirname(opt$outputpre)
  if(!file.exists(outdir)) dir.create(outdir)
}

pcasl <- tryCatch({
    antsImageRead(as.character(opt$pCASL), 4)
  }, error = function(e) {
    stop(paste('pCASL image', as.character(opt$pCASL),
                'does not exist.'))
})

if(length(opt$paramFile) > 0){
  if(file.exists(as.character(opt$paramFile))) {
    config <- read.csv(opt$paramFile)
  } else {
    config <- data.frame(tagFirst=T, sequence='pcasl')
  }
}

if (opt$smoothingFWHM > 0) {
  mysmoother <- c(rep(opt$smoothingFWHM, 3), 0)
  pcasl <- smoothImage(pcasl, mysmoother, FWHM=TRUE)
}

avg <- getAverageOfTimeSeries(pcasl)
avg <- n3BiasFieldCorrection(avg, 2)
avg <- n3BiasFieldCorrection(avg, 2)
mask <- getMask(avg, mean(avg), Inf, 2)
avg[mask==0] <- 0

moco <- antsrMotionCalculation(pcasl, fixed=avg, mask=mask)
tag.first <- config$tagFirst
ts <- timeseries2matrix(moco$moco_img, moco$moco_mask)
if (!tag.first) {
  tc <- (rep(c(1, 0), dim(ts)[1])[1:dim(ts)[1]] - 0.5)  # control minus tag
} else {
  tc <- (rep(c(0, 1), dim(ts)[1])[1:dim(ts)[1]] - 0.5)  # tag minus control
}
nuisance <- getASLNoisePredictors(ts, tc, polydegree='loess')
noise.all <- cbind(moco$moco_params, moco$fd$MeanDisplacement, nuisance)
noise.combined <- as.matrix(combineNuisancePredictors(ts, tc, noise.all))
onlypairs <- FALSE
if (opt$method == 'subtract') {
  onlypairs <- TRUE
}
censored <- aslCensoring(pcasl, mask, nuis=noise.combined, method='robust',
                         reject.pairs=onlypairs)
if (length(censored$which.outliers) > 0) {
  tc <- tc[-censored$which.outliers]
  noise.censored <- noise.combined[-censored$which.outliers, ]
} else {
  noise.censored <- noise.combined
}

if (opt$debug) {
  mean.ts <- apply(ts, 1, mean)
  dat.debug <- cbind(data.frame(MeanTimeSeries=mean.ts), noise.all)
  write.csv(dat.debug, file=paste(opt$outputpre, 'TimeSeriesData.csv', sep=''),
    row.names=as.character(1:nrow(ts)))
  write.csv(data.frame(Outliers=censored$which.outliers),
            file=paste(opt$outputpre, 'OutlierTimepoints.csv', sep=''))
}

if (opt$method == 'regression') {
  perf <- aslAveraging(censored$asl.inlier, mask=moco$moco_mask,
                tc=tc, nuisance=noise.censored, method='regression')
} else if (opt$method == 'bayesian') {
  if (length(opt$antsCorticalThicknessPrefix) == 0) {
    stop("For Bayesian regression, segmentations are required.")
  }
  act <- as.character(opt$antsCorticalThicknessPrefix)
  braint1 <- tryCatch({
        antsImageRead(paste(act, "ExtractedBrain0N4.nii.gz", sep=""))
      }, error = function(e) {
        print(paste('T1 brain image', paste(act, "ExtractedBrain0N4.nii.gz", sep=""),
                    'does not exist.'))
  })
  segmentation <- tryCatch({
        antsImageRead(paste(act, "BrainSegmentation.nii.gz", sep=""))
      }, error = function(e) {
        stop(paste('Segmentation image', paste(act, "BrainSegmentation.nii.gz", sep=""),
                    'does not exist.'))
  })
  postnames <- list.files(path=dirname(act),
    glob2rx("*BrainSegmentationPosteriors*.nii.gz"), full.names=TRUE)
  tissuelist <- tryCatch({
    imageFileNames2ImageList(postnames)
  }, error = function(e) {
    stop(paste("Probability images", postnames, "cannot be loaded."))
  })
  reg.t12asl <- antsRegistration(fixed=avg, moving=braint1,
    typeofTransform="SyNBold", outprefix=as.character(opt$outputpre))
  seg.asl <- antsApplyTransforms(avg, segmentation, reg.t12asl$fwdtransforms,
                                 "MultiLabel")
  for (ii in 1:length(tissuelist)) {
    tissuelist[[ii]] <- antsApplyTransforms(avg, tissuelist[[ii]],
                                    reg.t12asl$fwdtransforms, "Linear")
  }
  perf <- aslAveraging(censored$asl.inlier, mask=moco$moco_mask,
                tc=tc, nuisance=noise.censored, method='bayesian',
                       segmentation=seg.asl, tissuelist=tissuelist)
} else if(opt$method == 'subtract'){
  perf <- aslAveraging(censored$asl.inlier, mask=moco$moco_mask,
                       tc=tc, method='cubicSubtract')
}

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
m0 <- n3BiasFieldCorrection(m0,4)
m0 <- n3BiasFieldCorrection(m0,2)

if (length(opt$config > 0)) {
  tryCatch({
    config <- read.csv(opt$config, row.names=1)
  }, error = function(e){
    print(paste("Configuration file", opt$config, "does not exist."))
  })
  parameters <- c(list(m0=antsImageClone(m0)), config)
} else {
  parameters = list(sequence="pcasl", m0=antsImageClone(m0))
}

if (opt$debug) {
  antsImageWrite(perf, paste(opt$outputpre, 'Perfusion.nii.gz', sep=''))
  antsImageWrite(m0, paste(opt$outputpre, 'M0.nii.gz', sep=''))
}
cbf <- quantifyCBF(perf, mask=moco$moco_mask, parameters=parameters)
antsImageWrite(cbf$meancbf, paste(opt$outputpre, "CBF.nii.gz", sep=""))

if (nchar(opt$antsCorticalThicknessPrefix) > 0){
  act <- as.character(opt$antsCorticalThicknessPrefix)
  braint1 <- tryCatch({
      antsImageRead(paste(act, "ExtractedBrain0N4.nii.gz", sep=""))
    }, error = function(e) {
      print(paste('T1 brain image', paste(act, "ExtractedBrain0N4.nii.gz", sep=""),
                  'does not exist.'))
  })
  seg <- tryCatch({
      antsImageRead(paste(act, "BrainSegmentation.nii.gz", sep=""))
    }, error = function(e) {
      print(paste('Segmentation image', paste(act, "BrainSegmentation.nii.gz", sep=""),
                  'does not exist.'))
  })
  reg.t12asl <- antsRegistration(fixed=avg, moving=braint1,
    typeofTransform="SyNBold" )
  seg.asl <- antsApplyTransforms(avg, seg, reg.t12asl$fwdtransforms, "MultiLabel")
  antsImageWrite(seg.asl, paste(opt$outputpre,
                      "SegmentationWarpedToASL.nii.gz", sep=''))
  segstats <- labelStats(cbf$meancbf, seg.asl)
  write.csv(segstats, paste(opt$outputpre, 'TissueStats.csv', sep=''),
            row.names=FALSE)

  tx.template2t1 <- c(paste(act, "TemplateToSubject0Warp.nii.gz", sep=""),
                      paste(act, "TemplateToSubject1GenericAffine.mat", sep=""))
  tx.t12template <- c(paste(act, "SubjectToTemplate1Warp.nii.gz", sep=""),
                      paste(act, "SubjectToTemplate0GenericAffine.mat", sep=""))
  tx.asl2template <- c(reg.t12asl$invtransforms, tx.t12template)

  if (length(opt$template) > 0) {
    template <- tryCatch({
      antsImageRead(as.character(opt$template))
    }, error = function(e) {
      print(paste("Template image", template, "does not exist."))
    })
    asl.warped2template <- antsApplyTransforms(template, cbf$meancbf, tx.asl2template,
                                    whichtoinvert=c(F, F, F, F))
    antsImageWrite(asl.warped2template,
                   paste(opt$outputpre, "CBFWarpedToTemplate.nii.gz", sep=''))
  }

  tx.template2asl <- c(tx.template2t1, reg.t12asl$fwdtransforms)
  if (nchar(as.character(opt$labelSet)) > 0) {
    label <- tryCatch({
      antsImageRead(as.character(opt$labelSet))
    }, error = function(e) {
      print(paste("Label image", opt$labelSet, "does not exist."))
    })
    label.asl <- antsApplyTransforms(avg, label, tx.template2asl, "MultiLabel")
    antsImageWrite(label.asl, paste(opt$outputpre,
      'LabelWarpedToASL.nii.gz', sep=''))
    labelstats.cbf <- labelStats(cbf$meancbf, label.asl)
    write.csv(labelstats.cbf, paste(opt$outputpre, 'LabelStats.csv', sep=''),
              row.names=FALSE)
  }
}
