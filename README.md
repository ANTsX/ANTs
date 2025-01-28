[![ci-docker](https://github.com/ANTsX/ANTs/actions/workflows/ci-docker.yml/badge.svg)](https://github.com/ANTsX/ANTs/actions/workflows/ci-docker.yml)
[![Docker Pulls](https://img.shields.io/docker/pulls/antsx/ants.svg)](https://hub.docker.com/repository/docker/antsx/ants)
![Downloads](https://img.shields.io/github/downloads/antsx/ants/total)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/ants/badges/version.svg)](https://anaconda.org/conda-forge/ants)
[![PubMed](https://img.shields.io/badge/ANTsX_paper-Open_Access-8DABFF?logo=pubmed)](https://pubmed.ncbi.nlm.nih.gov/33907199/)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](code_of_conduct.md)

![ants template](http://i.imgur.com/mLZ71Ai.png)

**Advanced Normalization Tools (ANTs)** is a C++ library available through the command line that computes high-dimensional mappings to capture the statistics of brain structure and function. It allows one to organize, visualize and statistically explore large biomedical image sets. Additionally, it integrates imaging modalities in space + time and works across species or organ systems with minimal customization. 

The ANTs library is considered a state-of-the-art medical image registration and segmentation toolkit which depends on the Insight ToolKit, a widely used medical image processing library to which ANTs developers contribute. ANTs-related tools have also won several international, unbiased competitions such as MICCAI, BRATS, and STACOM.

It is possible to use ANTs in R ([ANTsR](https://github.com/antsx/antsr)) and Python ([ANTsPy](https://github.com/antsx/antspy)), with additional functionality for deep learning in R ([ANTsRNet](https://github.com/antsx/antsrnet)) and Python ([ANTsPyNet](https://github.com/antsx/antspynet)). These libraries help integrate ANTs with the broader R / Python ecosystem.

<br />

## Installation

Quick links: [download binaries](https://github.com/ANTsX/ANTs/releases) | [build from source](https://github.com/ANTsX/ANTs/wiki/Compiling-ANTs-on-Linux-and-Mac-OS) | [docker](https://hub.docker.com/r/antsx/ants) | [conda]((https://anaconda.org/conda-forge/ants).

### Pre-compiled binaries

The easiest way to install ANTs is by downloading the latest binaries on the [Releases](https://github.com/ANTsX/ANTs/releases) page. Download the latest release under the "Assets" section, then unzip the archive. Next, add the ANTs library to your PATH:

```
export PATH=/path/to/ants/bin:$PATH
```

You can check that this worked by running a command to find the path to any ANTs function:

```
which antsRegistration
```


If that works, you should be able to use the full functionality of ANTs from the command line or bash. You may wish to control multi-threading by setting the environment variable `ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS`.

### Building from source

When necessary, you can also build ANTs from the latest source code. A minimal example on Linux / Mac looks like this:

```bash
workingDir=${PWD}
git clone https://github.com/ANTsX/ANTs.git
mkdir build install
cd build
cmake \
    -DCMAKE_INSTALL_PREFIX=${workingDir}/install \
    ../ANTs 2>&1 | tee cmake.log
make -j 4 2>&1 | tee build.log
cd ANTS-build
make install 2>&1 | tee install.log
```

More details and a full downloadable installation script can be found in the [Linux/MacOS Guide](https://github.com/ANTsX/ANTs/wiki/Compiling-ANTs-on-Linux-and-Mac-OS). Building from source will generally work on Windows as well with some additional steps explained in the [Windows Guide](https://github.com/ANTsX/ANTs/wiki/Compiling-ANTs-on-Windows-10). Alternatively, it is also possible to install ANTs via [Docker](https://hub.docker.com/r/antsx/ants) or [Conda](https://anaconda.org/conda-forge/ants).

<br />

## Code examples

ANTs is a flexible library that can be used for a variety of applications and areas. Below is a collection of example scripts that - with a little effort - can be adapted to fit your specific needs. Some examples also include code for ANTsR or ANTsPy.

### Registration

- Basic registration [[Link](https://github.com/stnava/ANTs/blob/master/Scripts/newAntsExample.sh)]
- Basic registration with mask [[Link](https://github.com/ntustison/antsRegistrationWithMaskExample)]
- Large deformation [[Link](http://stnava.github.io/C/)]
- Asymmetry [[Link](http://stnava.github.io/asymmetry/)]
- Automobile registration [[Link](http://stnava.github.io/cars/)]
- Point-set mapping [[Link](http://stnava.github.io/chicken/)]
- Global optimization [[Link](http://stnava.github.io/butterfly/)]
  
### Template construction

- Brain template [[Link](http://ntustison.github.io/TemplateBuildingExample/)]
- Single subject template [[Link](https://github.com/ntustison/SingleSubjectTemplateExample)]
- "Cooking" tissue priors for templates [[Link](https://github.com/ntustison/antsCookTemplatePriorsExample)]
  
### Cortical thickness

- Basic cortical thickness [[Link](https://github.com/ntustison/antsCorticalThicknessExample)]
- Chimpanzee example [[Link](https://github.com/stnava/WHopkinsNHP/)]

### Segmentation
- N4 bias correction + Atropos [[Link](https://github.com/ntustison/antsAtroposN4Example)]
- Brain tumor segmentation [[Link](https://github.com/ntustison/BRATS2013/tree/master/SimpleExample)]

### Brain

- Basic brain mapping [[Link](http://stnava.github.io/BasicBrainMapping/)]
- Brain extraction [[Link](https://github.com/ntustison/antsBrainExtractionExample)]
- Multi-atlas joint label/intensity fusion [[Link](https://github.com/ntustison/MalfLabelingExample), [Link](https://github.com/qureai/Multi-Atlas-Segmentation)] (credit: @chsasank)
- fMRI or Motion Correction [[Link](http://stnava.github.io/fMRIANTs/)]
- fMRI reproducibility [[Link](http://stnava.github.io/RfMRI/)]
- Partial EPI slab to T1 image registration [[Link](https://github.com/ntustison/PartialSlabEpiT1ImageRegistration)]
  
See also our pre-built ANTs templates with spatial priors available for download [[General](http://figshare.com/articles/ANTs_ANTsR_Brain_Templates/915436), [MNI](https://figshare.com/articles/ANTs_files_for_mni_icbm152_nlin_sym_09a/8061914)].
  
### Lung

- CT lung registration [[Link](https://github.com/ntustison/antsCtLungRegistrationExample)]
- Lung mask registration [[Link](https://github.com/ntustison/ProtonCtLungMaskRegistration)]
- Lung and lobe estimation [[Link](https://github.com/ntustison/LungAndLobeEstimationExample)]
- Lung ventilation-based segmentation [[Link](https://github.com/ntustison/LungVentilationSegmentationExample)]

### Cardiac

- Basic example [[Link](http://stnava.github.io/LabelMyHeart)]
  
### Other

- Patch-based super-resolution [[Link](https://github.com/ntustison/NonLocalSuperResolutionExample)]
- Image denoising [[Link](https://github.com/ntustison/DenoiseImageExample)]
- Morphing [[Link](http://stnava.github.io/Morpheus/)]

<br />

## Learning resources

There are many different resources for learning about how to use ANTs functions and the methodology behind them. A selected list of useful resources is provided here.

* ANTs Wiki [[Link](https://github.com/ANTsX/ANTs/wiki)]
* ANTs Documentation [[Link](https://github.com/stnava/ANTsDoc/blob/master/ants2.pdf)]
* ANTs Tutorials [[Link](https://github.com/stnava/ANTsTutorial)]

Some commonly visited tutorials for specific ANTs functions are also presented below.

* Using antsRegistration [[Link](https://github.com/ANTsX/ANTs/wiki/ANTS-and-antsRegistration)]
* Applying warps with antsApplyTransforms [[Link](https://github.com/ANTsX/ANTs/wiki/Forward-and-inverse-warps-for-warping-images,-pointsets-and-Jacobians)]
* Using antsCorticalThickness [[Link](https://github.com/ANTsX/ANTs/wiki/antsCorticalThickness-and-antsLongitudinalCorticalThickness-output)]
* Using N4BiasFieldCorrection [[Link](https://github.com/ANTsX/ANTs/wiki/N4BiasFieldCorrection)]
* Multi-modality Presentation [[Link](https://github.com/stnava/ANTS_MultiModality/blob/master/ants_multimodality.pdf)]
  
<br />

## Contributing

If you have a question, feature request, or bug report the best way to get help is by posting an issue on the GitHub page. Please remember that it is difficult to provide any help if you do not provide enough information to reproduce your issue or environment.

We welcome any new contributions and ideas to improve ANTs. If you want to contribute code, the best way to get started is by reading through the [Wiki](https://github.com/ANTsX/ANTs/wiki) to get an understanding of the project or by posting an issue.

<br />

## Team

Development of ANTs is led by [Brian B. Avants](http://stnava.github.io/Resume/) (Creator, Algorithm Design, Implementation), [Nicholas J. Tustison](http://ntustison.github.io/CV/) (Compeller, Algorithm Design, Implementation Guru), Hans J. Johnson (Large-Scale Application, Testing, Software design), Gang Song (Originator), Philip A. Cook, Jeffrey T. Duda (DTI), Ben M. Kandel (Perfusion, multivariate analysis), and Nick Cullen (Python, R). 

<br />

## References

A large collection of journal articles have been published using ANTs software and can be found by searching Google Scholar or PubMed. Below, we provide a curated list of the most relevant articles to be used as a guide for better understanding or citing ANTs.

### Image Registration

<i>Symmetric diffeomorphic image registration with cross-correlation: evaluating automated labeling of elderly and neurodegenerative brain</i>.  Med Image Anal (2008). [[Link](http://www.ncbi.nlm.nih.gov/pubmed/17659998)]

<i>Evaluation of 14 nonlinear deformation algorithms applied to human brain MRI registration</i>. Neuroimage (2009). [[Link](http://www.ncbi.nlm.nih.gov/pubmed/19195496)]

<i>Evaluation of registration methods on thoracic CT: the EMPIRE10 challenge</i>. IEEE Trans Med Imaging (2011). [[Link](http://www.ncbi.nlm.nih.gov/pubmed/21632295)]

<i>A reproducible evaluation of ANTs similarity metric performance in brain image registration</i>. Neuroimage (2011). [[Link](https://pubmed.ncbi.nlm.nih.gov/20851191/)]

### Templates

<i>The optimal template effect in hippocampus studies of diseased populations</i>. Neuroimage (2010). [[Link](https://pubmed.ncbi.nlm.nih.gov/19818860/)]

### Image Segmentation

<i>An open source multivariate framework for n-tissue segmentation with evaluation on public data</i>. Neuroinformatics (2011). [[Link](http://www.ncbi.nlm.nih.gov/pubmed/21373993)]

<i>Multi-atlas segmentation with joint label fusion and corrective learning—an open source implementation</i>. Front Neuroinform (2013). [[Link](https://www.frontiersin.org/articles/10.3389/fninf.2013.00027/full)]

### Bias Correction

<i>N4ITK: improved N3 bias correction</i>. IEEE Trans Med Imaging (2010). [[Link](http://www.ncbi.nlm.nih.gov/pubmed/20378467)]

### Cortical Thickness

<i>Registration based cortical thickness measurement</i>. Neuroimage (2009). [[Link](http://www.ncbi.nlm.nih.gov/pubmed/19150502)]

<i>Large-scale evaluation of ANTs and FreeSurfer cortical thickness measurements</i>. Neuroimage (2014). [[Link](https://pubmed.ncbi.nlm.nih.gov/24879923/)]

<i>Regional and hemispheric variation in cortical thickness in chimpanzees</i>. J Neurosci (2013). [[Link](http://www.ncbi.nlm.nih.gov/pubmed/23516289)]

<i>Longitudinal Mapping of Cortical Thickness Measurements: An Alzheimer's Disease Neuroimaging Initiative-Based Evaluation Study</i>.  J Alzheimers Dis
 (2019). [[Link](https://pubmed.ncbi.nlm.nih.gov/31356207/)]
 
### Eigenanatomy 

<i>Eigenanatomy improves detection power for longitudinal cortical change</i>. Med Image Comput Comput Assist Interv (2012). [[Link](http://www.ncbi.nlm.nih.gov/pubmed/23286132)]

<i>White matter imaging helps dissociate tau from TDP-43 in frontotemporal lobar degeneration</i>. J Neurol Neurosurg Psychiatry (2013). [[Link](https://pubmed.ncbi.nlm.nih.gov/23475817/)]

### Software

<i>The ANTsX ecosystem for quantitative biological and medical imaging</i>. Scientific Reports (2021). [[Link](https://www.nature.com/articles/s41598-021-87564-6)]

<i>ANTsX neuroimaging-derived structural phenotypes of UK Biobank</i>.  Scientific Reports (2024). [[Link](https://pubmed.ncbi.nlm.nih.gov/38632390/)]
<br />

## Funding

Current support comes from R01-EB031722. Previous support includes R01-EB006266-01 and K01-ES025432-01.
