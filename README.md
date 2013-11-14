============================
Advanced Normalization Tools 
============================

Questions: [Discussion Site](http://sourceforge.net/p/advants/discussion/)

 Email: [antsr.me at gmail dot com](mailto:antsr.me@gmail.com)

 Install ANTs via pre-built:
[Packages](http://sourceforge.net/projects/advants/files/ANTS/)

 Build [ANTs](https://github.com/stnava/ANTs) from:
[Source-Code](http://brianavants.wordpress.com/2012/04/13/updated-ants-compile-instructions-april-12-2012/)
(recommended)

ANTs extracts information from complex datasets that include imaging
([Word Cloud](http://brianavants.files.wordpress.com/2013/05/avants_wordcloud.jpg)).
Paired with [ANTsR](http://stnava.github.com/ANTsR/) (answer), ANTs is
useful for managing, interpreting and visualizing multidimensional data.
ANTs is
[popularly](https://sourceforge.net/projects/advants/files/ANTS/stats/timeline?dates=2010-07-19+to+2099-05-25)
considered a state-of-the-art medical image registration and
segmentation toolkit. ANTsR is an emerging tool supporting standardized
multimodality image analysis. ANTs depends on the Insight ToolKit
[(ITK)](http://www.itk.org), a widely used medical image processing
library to which ANTs developers contribute.

Authors
-------

### Brian B. Avants - UPENN

**Role:** Creator, Algorithm Design, Implementation,
[more](http://stnava.github.io/Resume/)

### Nicholas J. Tustison - UVA

**Role:** Compeller, Algorithm Design, Implementation Guru

### Hans J. Johnson - UIowa

**Role:** Large-Scale Application, Testing, Software design

### Team Members 

**Core:** Gang Song (Originator), Jeffrey T. Duda (DTI), Ben M. Kandel (Perfusion, multivariate analysis), Kent Williams (software engineer, UIowa) 


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
methods [(link)](http://www.ncbi.nlm.nih.gov/pubmed/21237273N4), Bias
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

-   New *antsRegistration* [bash
    example](https://github.com/stnava/ANTs/blob/master/Scripts/newAntsExample.sh)

-   New *Eigenanatomy* for [multivariate neuroimage
    analysis](http://www.ncbi</a>.nlm.nih.gov/pubmed/23269595) via
    [PCA](http://www.ncbi.nlm.nih.gov/pubmed/23286132) &
    [CCA](http://www.ncbi.nlm.nih.gov/pubmed/20083207)

-   ANTs finished in 1st rank in [Klein 2009 intl. brain mapping
    competition](http://www.ncbi.nlm.nih.gov/pubmed/19195496)

-   ANTs finished 1st overall in [EMPIRE10 intl. lung mapping
    competition](http://www.ncbi.nlm.nih.gov/pubmed/21632295)

-   ANTs is the standard registration for
    [MICCAI-2013](http://www.miccai2013.org/) segmentation competitions

-   Conducting ANTs-based R tutorial @ MICCAI-2013

-   ITK-focused Frontiers in Neuroinformatics research topic
    [here](http://www.frontiersin.org/neuroinformatics/researchtopics/neuroinformatics_with_the_insi/1580)

-   Nick led us to a win at the [BRATS 2013 challenge](http://martinos.org/qtim/miccai2013/) with [ANTsR](http://stnava.github.io/ANTsR/)

Learning about ANTs
----------------------------------------------------------

**Basic Brain Mapping** [example](http://stnava.github.io/BasicBrainMapping/)

**Large deformation** [example](http://stnava.github.io/C/)

**Template construction** [example](https://github.com/ntustison/TemplateBuildingExample)

**Automobile** [example](http://stnava.github.io/cars/)

**Asymmetry** [example](http://stnava.github.io/asymmetry/)

**Point-set** [mapping](http://stnava.github.io/chicken/)

**Feature matching** [example](http://stnava.github.io/featureMatching/) ... not up to date ...

**Chimpanzee cortical thickness** [example](http://stnava.github.io/WHopkinsNHP/)

**Global optimization** [example](http://stnava.github.io/butterfly/)

**Morphing** [example](http://stnava.github.io/Morpheus/)

**fMRI or Motion Correction** [example](http://stnava.github.io/fMRIANTs/)

**fMRI reproducibility** [example](http://stnava.github.io/RfMRI/)

**fMRI prediction** [example](http://stnava.github.io/Haxby2001/) ... WIP ...

**Cardiac** [example](http://stnava.github.io/LabelMyHeart)

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

References
----------------------------------------------------------

[Google
Scholar](http://scholar.google.com/scholar?q=Advanced+Normalization+Tools+%22ANTs%22+-ant&hl=en&as_sdt=1%2C39&as_ylo=2008&as_yhi=)

[Pubmed](http://www.ncbi.nlm.nih.gov/pubmed?term=%22Avants%20B%22%20OR%20%22Tustison%20N%22)

ANTs was supported by: R01-EB006266-01
