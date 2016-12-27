<!--
![ants registration artillery](http://i.imgur.com/FCLrXV1.jpg)
![ants multivar](http://i.imgur.com/YqWtunL.png)
![ants faces](http://i.imgur.com/wBOFGwg.png)
![ants goat](http://i.imgur.com/SEKf1mo.jpg) -->
============================
[Advanced Normalization Tools](https://imgur.com/a/kySGi)
============================

[![Build Status](https://travis-ci.org/stnava/ANTs.svg?branch=master)](https://travis-ci.org/stnava/ANTs)

ANTs computes high-dimensional mappings to capture the statistics of brain
structure and function.  See the [FAQ page](https://github.com/stnava/ANTsTutorial/blob/master/handout/antsGithubExamples.Rmd).

![ants template](http://i.imgur.com/mLZ71Ai.png)

ANTs allows one to organize, visualize and statistically explore large biomedical
image sets.

![ants render](http://i.imgur.com/hMW6fjB.png)

ANTs integrates imaging modalities and related information in space and time.

![ants render](http://i.imgur.com/oIMrnpY.png)

ANTs works across species or organ systems with minimal customization.

![ants primate](http://i.imgur.com/Dfrifgg.png)

ANTs and related tools have won several international and unbiased competitions.

![ants competes](http://i.imgur.com/HE0j7IC.png)

[ANTsR](http://stnava.github.io/ANTsR/) is the underlying statistical workhorse.

Questions: [Discussion Site](http://sourceforge.net/p/advants/discussion/) or *new* [ANTsDoc](http://stnava.github.io/ANTsDoc/) or try [this version](http://issuu.com/brianavants/docs/ants2) ... also read our [guide to evaluation strategies and addressing new problems with ANTs or other software](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3766821/).  *New* [ANTs handout](https://github.com/stnava/ANTsTutorial/raw/master/handout/antsHandout.pdf), part of forthcoming [ANTs tutorial]() material.

[ANTsTalk - subject to change at any moment](http://stnava.github.io/ANTsTalk/)

[ANTsRegistrationTalk - subject to change at any moment](http://stnava.github.io/ANTsRegistrationTalk/)

 Install ANTs via pre-built:
[Packages @ github](https://github.com/stnava/ANTs/releases) older
versions [@ sourceforge](http://sourceforge.net/projects/advants/files/ANTS/) ... also,
[Github Releases are here](https://github.com/stnava/ANTs/releases) thanks to Arman Eshaghi.

 Build [ANTs](https://github.com/stnava/ANTs) from:
[Source-Code](http://brianavants.wordpress.com/2012/04/13/updated-ants-compile-instructions-april-12-2012/)
(recommended)

[ANTs Dashboard](https://travis-ci.org/stnava/ANTs/) thanks to Arman Eshaghi and  Hans J. Johnson

ANTs extracts information from complex datasets that include imaging
([Word Cloud](http://brianavants.files.wordpress.com/2013/05/avants_wordcloud.jpg)).
Paired with [ANTsR](http://stnava.github.io/software/2014/01/08/antsr/) (answer), ANTs is
useful for managing, interpreting and visualizing multidimensional data.
ANTs is
[popularly](https://sourceforge.net/projects/advants/files/ANTS/stats/timeline?dates=2010-07-19+to+2099-05-25)
considered a state-of-the-art medical image registration and
segmentation toolkit. ANTsR is an emerging tool supporting standardized
multimodality image analysis. ANTs depends on the Insight ToolKit
[(ITK)](http://www.itk.org), a widely used medical image processing
library to which ANTs developers contribute.  A summary of some ANTs findings and tutorial material (most of which is on this page) is [here](http://rpubs.com/stnava/ANTsTut).

Authors
-------

### Brian B. Avants - UPENN

**Role:** Creator, Algorithm Design, Implementation,
[more](http://stnava.github.io/Resume/)

### Nicholas J. Tustison - UVA

**Role:** Compeller, Algorithm Design, Implementation Guru,
[more](http://ntustison.github.io/CV/)

### Hans J. Johnson - UIowa

**Role:** Large-Scale Application, Testing, Software design

### Team Members

**Core:** Gang Song (Originator), Philip A. Cook, Jeffrey T. Duda (DTI), Ben M. Kandel (Perfusion, multivariate analysis)

Image Registration
------------------

Diffeomorphisms: [SyN](http://www.ncbi.nlm.nih.gov/pubmed/17659998),
Independent Evaluation:
[Klein](http://www.ncbi.nlm.nih.gov/pubmed/19195496),
[Murphy](http://www.ncbi.nlm.nih.gov/pubmed/21632295), Template
Construction
[(2004)](http://www.ncbi.nlm.nih.gov/pubmed/15501083)[(2010)](http://www.ncbi.nlm.nih.gov/pubmed/19818860),
[Similarity Metrics](http://www.ncbi.nlm.nih.gov/pubmed/20851191),
[Multivariate
registration](http://www.ncbi.nlm.nih.gov/pubmed/18995188), [Multiple
modality analysis and statistical
bias](http://www.ncbi.nlm.nih.gov/pubmed/23151955)

Image Segmentation
------------------

Atropos Multivar-EM Segmentation
[(link)](http://www.ncbi.nlm.nih.gov/pubmed/21373993), Multi-atlas
methods [(link)](https://scholar.google.com/scholar?q=joint+label+fusion+yushkevich&btnG=&hl=en&as_sdt=0%2C31) and [JLF](http://journal.frontiersin.org/article/10.3389/fninf.2013.00027/full), Bias
Correction [(link)](http://www.ncbi.nlm.nih.gov/pubmed/20378467), DiReCT
cortical thickness
[(link)](http://www.ncbi.nlm.nih.gov/pubmed/19150502), DiReCT in
[chimpanzees](http://www.ncbi.nlm.nih.gov/pubmed/23516289)

Multivariate Analysis Eigenanatomy [(1)](http://www.ncbi.nlm.nih.gov/pubmed/23286132) [(2)](http://www.ncbi.nlm.nih.gov/pubmed/23475817)
----------------------------------------------------------------------------------------------------------------------------------------

Prior-Based Eigenanatomy [(in
prep)](http://www.ncbi.nlm.nih.gov/pubmed/?), Sparse CCA
[(1)](http://www.ncbi.nlm.nih.gov/pubmed/20083207),
[(2)](http://www.ncbi.nlm.nih.gov/pubmed/20879247), Sparse Regression
[(link)](http://link.springer.com/chapter/10.1007%2F978-3-642-38868-2_8)

ImageMath Useful!
-----------------

morphology, GetLargestComponent, CCA, FillHoles ... much more!

Application Domains
-------------------

### Frontotemporal degeneration [PENN FTD center](http://ftd.med.upenn.edu)

### Multimodality Neuroimaging

-   [Structural MRI](http://jeffduda.github.io/NeuroBattery/)
-   Functional MRI
-   Network Analysis

### Lung Imaging

-   Structure
-   Perfusion MRI
-   Branching

### Multiple sclerosis (lesion filling) [example](https://github.com/armaneshaghi/LesionFilling_example)

Background & Theory
----------------------------------------------------------

-   The
    [SyN](http://www.ncbi.nlm.nih.gov/pubmed/?term=%22SyN%22+AND+%22Avants+B%22)
    and [N4 bias
    correction](http://www.ncbi.nlm.nih.gov/pubmed/?term=%22N4%22+AND+%22Tustison+N4ITK%22)
    papers and other relevant references in
    [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/?term=%22Tustison+N%22+AND+%22Avants+B%22)

-   Visualization: e.g. [a gource of ANTs
    development](http://vimeo.com/66781467)

-   [DiReCT](http://www.ncbi.nlm.nih.gov/pubmed/?term=%22DIRECT%22+AND+%22Avants%22+AND+DAS)
    cortical thickness
    [papers](http://www.ncbi.nlm.nih.gov/pubmed/?term=%22Cortical+Thickness%22+AND+%22Avants%22)

-   A
    [folder](https://sourceforge.net/projects/advants/files/Documentation/)
    of relevant docs:
    [segmentation](http://sourceforge.net/projects/advants/files/Documentation/atropos.pdf/download),
    [registration](http://sourceforge.net/projects/advants/files/Documentation/antstheory.pdf/download),
    [usage(old)](http://sourceforge.net/projects/advants/files/Documentation/ants.pdf/download),
    [for clinical
    apps](http://sourceforge.net/projects/advants/files/Documentation/ANTSMethodologySummary.docx/download)

-   ANTs redesigned for generality, automation, multi-core computation
    with ITKv4

-   Dev'd ITKv4 with Kitware, GE, Natl. Lib of Medicine & Academia


ANTs has won several unbiased & international competitions
----------------------------------------------------------

-   ANTs finished in 1st rank in [Klein 2009 intl. brain mapping
    competition](http://www.ncbi.nlm.nih.gov/pubmed/19195496)

-   ANTs finished 1st overall in [EMPIRE10 intl. lung mapping
    competition](http://www.ncbi.nlm.nih.gov/pubmed/21632295)

-   ANTs is the standard registration for
    [MICCAI-2013](http://www.miccai2013.org/) segmentation competitions

-   Conducting ANTs-based R tutorial @ MICCAI-2013

-   ITK-focused Frontiers in Neuroinformatics research topic
    [here](http://www.frontiersin.org/neuroinformatics/researchtopics/neuroinformatics_with_the_insi/1580)

-   Won the [BRATS 2013 challenge](http://martinos.org/qtim/miccai2013/) with [ANTsR](http://stnava.github.io/ANTsR/)

-   Won the best paper award at the [STACOM 2014 challenge](http://www.cardiacatlas.org/web/stacom2014/home)

Learning about ANTs
----------------------------------------------------------
**antsRegistration** [bash example](https://github.com/stnava/ANTs/blob/master/Scripts/newAntsExample.sh)

**Eigenanatomy** for [multivariate neuroimage analysis](http://www.ncbi</a>.nlm.nih.gov/pubmed/23269595) via
    [PCA](http://www.ncbi.nlm.nih.gov/pubmed/23286132) &
    [CCA](http://www.ncbi.nlm.nih.gov/pubmed/20083207)

**ANTs and ITK** [paper](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4009425/)

**Pre-built ANTs templates with spatial priors** [download](http://figshare.com/articles/ANTs_ANTsR_Brain_Templates/915436)

**The ANTs Cortical Thickness Pipeline** [example](https://github.com/ntustison/KapowskiChronicles/blob/master/paper2.pdf?raw=true)

**"Cooking" tissue priors for templates**
  [example](https://github.com/ntustison/antsCookTemplatePriorsExample)
  (after you build your template)

**Basic Brain Mapping** [example](http://stnava.github.io/BasicBrainMapping/)

**Large deformation** [example](http://stnava.github.io/C/)

**Template construction** [example](http://ntustison.github.io/TemplateBuildingExample/)

**Automobile** [example](http://stnava.github.io/cars/)

**Asymmetry** [example](http://stnava.github.io/asymmetry/)

**Point-set** [mapping](http://stnava.github.io/chicken/) which includes the
PSE metric and affine and deformable registration with (labeled) pointsets or
iterative closest point

**Feature matching** [example](http://stnava.github.io/featureMatching/) ... not up to date ...

**Chimpanzee cortical thickness** [example](https://github.com/stnava/WHopkinsNHP/)

**Global optimization** [example](http://stnava.github.io/butterfly/)

**Morphing** [example](http://stnava.github.io/Morpheus/)

**fMRI or Motion Correction** [example](http://stnava.github.io/fMRIANTs/)

**fMRI reproducibility** [example](http://stnava.github.io/RfMRI/)

**fMRI prediction** [example](http://stnava.github.io/Haxby2001/) ... WIP ...

**Cardiac** [example](http://stnava.github.io/LabelMyHeart)

**Brain extraction** [example](https://github.com/ntustison/antsBrainExtractionExample)

**N4 bias correction <-> segmentation** [example](https://github.com/ntustison/antsAtroposN4Example)

**Cortical thickness** [example](https://github.com/ntustison/antsCorticalThicknessExample)

**Multi-atlas joint label/intensity fusion examples** [example 1](https://github.com/ntustison/MalfLabelingExample) [example 2](https://github.com/qureai/Multi-Atlas-Segmentation) (thanks to @chsasank)

**Lung and lobe estimation** [example](https://github.com/ntustison/LungAndLobeEstimationExample)

**Bibliography** [bibtex of ANTs-related papers](https://github.com/stnava/ANTsBibliography)

**ANTs** [google scholar page](http://scholar.google.com/citations?user=ox-mhOkAAAAJ&hl=en)

Presentations: e.g. [a Prezi about
ANTs](http://prezi.com/mwrmcm-h9-w4/ants/?kw=view-mwrmcm-h9-w4&rc=ref-40024395)
(WIP)

Reproducible science as a teaching tool: e.g. [compilable ANTs
tutorial](https://github.com/stnava/ANTS_MultiModality) (WIP)

Other examples [slideshow](http://brianavants.wordpress.com)

Landmark-based mapping for e.g. hippocampus [discussed
here](https://sourceforge.net/p/advants/discussion/840261/thread/1cb7b165/?limit=50)

Brief ANTs segmentation [video](http://vimeo.com/67814201)

**Benchmarks** for expected memory and computation time: [results](https://github.com/gdevenyi/antsRegistration-benchmarking).  These
results are, of course, system and data dependent.

References
----------------------------------------------------------

[Google
Scholar](http://scholar.google.com/scholar?q=Advanced+Normalization+Tools+%22ANTs%22+-ant&hl=en&as_sdt=1%2C39&as_ylo=2008&as_yhi=)

[Pubmed](http://www.ncbi.nlm.nih.gov/pubmed?term=%22Avants%20B%22%20OR%20%22Tustison%20N%22)

Boilerplate ANTs
------------------

Here is some boilerplate regarding ants image processing:

We will analyze multiple modality neuroimaging data with Advanced
Normalization Tools (ANTs) version >= 2.1 [1]
(http://stnava.github.io/ANTs/).  ANTs has proven performance in
lifespan analyses of brain morphology [1] and function [2] in both
adult [1] and pediatric brain data [2,5,6] including infants [7].
ANTs employs both probabilistic tissue segmentation (via Atropos [3])
and machine learning methods based on expert labeled data (via joint
label fusion [4]) in order to maximize reliability and consistency of
multiple modality image segmentation.  These methods allow detailed
extraction of critical image-based biomarkers such as volumes
(e.g. hippocampus and amygdala), cortical thickness and area and
connectivity metrics derived from structural white matter [13] or
functional connectivity [12]. Critically, all ANTs components are
capable of leveraging multivariate image features as well as expert
knowledge in order to learn the best segmentation strategy available
for each individual image [3,4].  This flexibility in segmentation and
the underlying high-performance normalization methods have been
validated by winning several internationally recognized medical image
processing challenges conducted within the premier conferences within
the field and published in several accompanying articles
[8][9][10][11].

References

[1] http://www.ncbi.nlm.nih.gov/pubmed/24879923

[2] http://www.ncbi.nlm.nih.gov/pubmed/24817849

[3] http://www.ncbi.nlm.nih.gov/pubmed/21373993

[4] http://www.ncbi.nlm.nih.gov/pubmed/21237273

[5] http://www.ncbi.nlm.nih.gov/pubmed/22517961

[6] http://www.ncbi.nlm.nih.gov/pubmed/24033570

[7] http://www.ncbi.nlm.nih.gov/pubmed/24139564

[8]  http://www.ncbi.nlm.nih.gov/pubmed/21632295

[9] http://www.ncbi.nlm.nih.gov/pubmed/19195496

[10] http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3837555/

[11] http://nmr.mgh.harvard.edu/~koen/MenzeTMI2014.pdf

[12] http://www.ncbi.nlm.nih.gov/pubmed/23813017

[13] http://www.ncbi.nlm.nih.gov/pubmed/24830834


ANTs was supported by: R01-EB006266-01 and by K01-ES025432-01

![ants chimp](http://i.imgur.com/4tPvy05.png)
