#!/usr/bin/env Rscript 
require(ANTsR)
usePkg('optparse')
usePkg('tools')

optlist <- list(make_option(c('-a', '--anatomical'), default='', help='anatomical image (skull stripped)'),
  make_option(c('-p', '--priors'), default='', help='brain segmentation priors (C-style, e.g. priors%d.nii.gz)'),
  make_option(c('-g', '--segmentation'), default='', help='hard brain segmentation'), 
  make_option(c('-x', '--mask'), default='', help='t1 brain mask'),
  make_option(c('-s', '--pCASL'), default='', help=' raw pCASL image'),
  make_option(c('-e', '--template'), default='', help='brain template'),
  make_option(c('-t', '--transformpre'), default='', help='skull-stripped t1 to template transform prefix'),
  make_option(c('-l', '--labels'),  default='', help='template labels'),
  make_option(c('-o', '--outputpre'), default='CBF_', help='output prefix (defaults to %default)'),
  make_option(c('-m', '--method'), default='RobustRegression', help=' method for perfusion calculation. \n\t\tOne of: "SimpleSubtraction", "SurroundSubtraction", "SincSubtraction", "RobustRegression", "BayesianRegression", "LocalBayesianRegression."'),
  make_option(c('-d', '--denoising'), default='CompCorMotion', help='denoising method.  Options are: "CompCor", "Motion", "Detrending", \n\t\t"Cross-Validation", "OutlierRejection".  Multiple options can be specified (e.g., "CompCorMotion" is legal).  Default is %default.'), 
  make_option(c('-b', '--bloodT1'), default=0.67, help='blood T1 value (defaults to %default s^-1)'),
  make_option(c('-r', '--robustness'), default=0.95, help='robustness parameter (defaults to %default).'),
  make_option(c('-n', '--bootstrapNumber'), default=20, help=' number of bootstrap samples (defaults to %default)'),
  make_option(c('-c', '--bootstrapPercent'), default=0.70, help='percent to sample per bootstrap run (defaults to %default)'),
  make_option(c('-k', '--keepTmp'), default=F, action='store_true', help='keep tmp files, including warps (defaults to %default--takes lots of space to save)'),
  make_option(c('-f', '--bootstrapReplace'), default=F, action='store_true', help='bootstrap with replacement?  takes arguments "false" or "true"; defaults to false.'), 
  make_option(c('-v', '--verbose'), default=F, action='store_true', help='verbose output.')) 

usage <- OptionParser(option_list=optlist, usage='Usage: %prog <apgxset> [lomdbrnckfv]') 
opt <- parse_args(usage) 

if(!file.exists(opt$anatomical)) 
  stop(paste('Anatomical image', opt$anatomical, 
    'does not exist.  For usage, see -h option.')
if(!file.exists(opt$segmentation))
  stop(paste('Brain segmentation', opt$segmentation, 
    'does not exist.  For usage, see -h option.')
if(!file.exists(opt$mask)) 
  stop(paste('Brain mask', opt$mask, 
    'does not exist.  For usage, see -h option.')
if(!file.exists(opt$template))
  stop(paste('Template', opt$template,  
    'does not exist.  For usage, see -h option.')
if(!file.exists(opt$pCASL)) 
  stop(paste('pCASL image', opt$pCASL, 
    'does not exist.  For usage, see -h option.')

if(opt$verbose) {
  cat('Running antsASLProcessing.R with the following options:\n')
  for(option in names(opt)){
    cat(paste(option, ': ', opt[option], '\n', sep=''))
  }
}

if(length(grep(.Platform$file.sep, opt$outputpre)) > 0) {
  outdir <- dirname(opt$outpre)
  if(!file.exists(outdir)) dir.create(outdir) 
}

problist <- list() 
for(ii in 1:6){
  probname <- sprintf(opt$priors, ii)
  if(!file.exists(probname)) 
    stop('6-tissue probability priors do not exist.  For usage, see -h option.')
  problist[[ii]] <- antsImageRead(probname, 3) 
}

template <- antsImageRead(opt$template, 3)
templateseg <- antsImageRead(opt$segmentation, 3) 
pcasl <- antsImageRead(opt$pCASL, 4)
avg <- getAverageOfTimeSeries(pcasl)
N3BiasFieldCorrection(3, avg, avg, 2)
N3BiasFieldCorrection(3, avg, avg, 2)
maskavg <- getMask(avg, mean(avg), Inf, 2)
avg[maskavg==0] <- 0

if(opt$verbose) print('Performing pCASL to template registration.')
reg.pcasl2template <- antsRegistration(avg, template, 
  typeofTransform='SyNBold',  outprefix=opt$outputpre)
seg.pcasl <- antsApplyTransforms(avg, templateseg,
   reg.pcasl2template$fwdtransforms,
   "MultiLabel") 
pcasl.probs <- list()
for (ii in 1:6){
  prob.warped2pcasl <- antsApplyTransforms(avg,
    problist[[ii]], reg.pcasl2template$fwdtransforms)
  pcasl.probs[[i]]<-prob.warped2pcasl
  antsImageWrite(prob.warped2pcasl, 
    paste(outfn, 'Probability0', i, '.nii.gz', sep=''))
}


