<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">

<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

<head>
     <meta http-equiv="Content-Type" content="text/html; charset=utf-8"/>

     <title>Advanced Normalization Tools</title>

     <style type="text/css">
        * { margin: 0; padding: 0; }
        body { font: 16px Helvetica, Sans-Serif; line-height: 24px; background: url(forhtml/nazca-1.jpg); }
        .clear { clear: both; }
        #page-wrap { width: 800px; margin: 40px auto 60px; }
        #pic { float: right; margin: -30px 0 0 0; }
        h1 { margin: 0 0 16px 0; padding: 0 0 16px 0; font-size: 42px; font-weight: bold; letter-spacing: -2px; border-bottom: 1px solid #999; }
        h2 { font-size: 20px; margin: 0 0 6px 0; position: relative; }
        h2 span { position: absolute; bottom: 0; right: 0; font-style: italic; font-family: Georgia, Serif; font-size: 16px; color: #00006E; font-weight: normal; }
        p { margin: 0 0 16px 0; }
        a { color: #980000 ; text-decoration: none; border-bottom: 1px dotted #00006E; }
        a:hover { border-bottom-style: solid; color: black; }
        ul { margin: 0 0 32px 17px; }
        #objective { width: 500px; float: left; }
        #objective p { font-family: Georgia, Serif; font-style: italic; color: #00006E; }
        dt { font-style: italic; font-weight: bold; font-size: 18px; text-align: right; padding: 0 26px 0 0; width: 150px; float: left; height: 100px; border-right: 1px solid #999;  }
        dd { width: 600px; float: right; }
        dd.clear { float: none; margin: 0; height: 15px; }
     </style>
</head>

<body>

    <div id="page-wrap">

        <div id="contact-info" class="vcard">

            <!-- Microformats! -->

# Advanced Normalization Tools

                Questions: [Discussion
                Site](http://sourceforge.net/p/advants/discussion/)

                Email: [antsr.me at gmail
                dot com](mailto:antsr.me@gmail.com)

	        Install ANTs via pre-built: [Packages](http://sourceforge.net/projects/advants/files/ANTS/)

	        Build [ANTs](https://github.com/stnava/ANTs) from: [Source-Code](http://brianavants.wordpress.com/2012/04/13/updated-ants-compile-instructions-april-12-2012/) (recommended)

        </div>

        <div id="objective">

            ANTs extracts information from
            complex datasets that include imaging ([Word
            Cloud](http://brianavants.files.wordpress.com/2013/05/avants_wordcloud.jpg)).  Paired with [ANTsR](http://stnava.github.com/ANTsR/) (answer),
            ANTs is useful for managing, interpreting and
            visualizing multidimensional data. ANTs is [popularly](https://sourceforge.net/projects/advants/files/ANTS/stats/timeline?dates=2010-07-19+to+2099-05-25) considered
            a state-of-the-art medical image registration and segmentation toolkit. ANTsR is an emerging tool supporting standardized
            multimodality image analysis.  ANTs depends on the Insight
            ToolKit  [(ITK)](http://www.itk.org), a widely used  
            medical image processing library to which ANTs developers contribute.

        </div>

        <div class="clear"></div>

        <dl>
            <dd class="clear"></dd>

            <dt>Authors</dt>
            <dd>

## Brian B. Avants - UPENN

**Role:** Creator, Algorithm Design,
            Implementation, [more](http://stnava.github.io/Resume/)

           </dd>
           <dd>

## Nicholas J. Tustison - UVA

**Role:** Compeller, Algorithm Design,
            Implementation Guru

            </dd>
	     <dd>

**Team:** Gang Song (Originator),
            Jeffrey T. Duda (DTI), Hans J. Johnson (Large-Scale
            Application, Testing)

            </dd>
            <dd class="clear"></dd>

            <dt>Methods References</dt>
            <dd>

## Image Registration

Diffeomorphisms: [SyN](http://www.ncbi.nlm.nih.gov/pubmed/17659998),
            Independent Evaluation: [Klein](http://www.ncbi.nlm.nih.gov/pubmed/19195496),
	      [Murphy](http://www.ncbi.nlm.nih.gov/pubmed/21632295),
	      Template Construction
	      [(2004)](http://www.ncbi.nlm.nih.gov/pubmed/15501083)[ (2010)](http://www.ncbi.nlm.nih.gov/pubmed/19818860),
	      [Similarity
            Metrics](http://www.ncbi.nlm.nih.gov/pubmed/20851191),
	      [Multivariate
            registration](http://www.ncbi.nlm.nih.gov/pubmed/18995188),
	      [Multiple
            modality analysis and statistical bias](http://www.ncbi.nlm.nih.gov/pubmed/23151955)

## Image Segmentation

Atropos Multivar-EM Segmentation [(link)](http://www.ncbi.nlm.nih.gov/pubmed/21373993),
            Multi-atlas methods [(link)](http://www.ncbi.nlm.nih.gov/pubmed/21237273N4), Bias
                Correction [(link)](http://www.ncbi.nlm.nih.gov/pubmed/20378467),
                DiReCT cortical thickness [(link)](http://www.ncbi.nlm.nih.gov/pubmed/19150502),
            DiReCT in
            [chimpanzees](http://www.ncbi.nlm.nih.gov/pubmed/23516289)

## Multivariate Analysis <span>Eigenanatomy [(1)](http://www.ncbi.nlm.nih.gov/pubmed/23286132)
                [(2)](http://www.ncbi.nlm.nih.gov/pubmed/23475817)
                </span>

Prior-Based Eigenanatomy [(in
                prep)](http://www.ncbi.nlm.nih.gov/pubmed/?),
                Sparse CCA [(1)](http://www.ncbi.nlm.nih.gov/pubmed/20083207),
	      [(2)](http://www.ncbi.nlm.nih.gov/pubmed/20879247),
                Sparse Regression [(link)](http://link.springer.com/chapter/10.1007%2F978-3-642-38868-2_8)

## ImageMath <span>Useful!</span>

morphology, GetLargestComponent, CCA, FillHoles
                ... much more!

          </dd>

            <dd class="clear"></dd>

            <dt>Application Domains</dt>
            <dd>

## Frontotemporal degeneration <span>[PENN
            FTD center](http://ftd.med.upenn.edu)</span>
               <h2>Neuroimaging <span>Multimodality</span>

*   Structural MRI
*   Functional MRI
*   Network Analysis 

## Lung imaging <span>Translational</span>

*   Structure
*   Perfusion MRI
*   Branching 
           </dd>

	    <dd class="clear"></dd>

            <dt>Background & Theory</dt>
            <dd>

*   The [
            SyN](http://www.ncbi.nlm.nih.gov/pubmed/?term=((%22SyN%22)+AND+%22Avants+B%22) and [
            N4 bias correction](http://www.ncbi.nlm.nih.gov/pubmed/?term=((%22N4%22)+AND+%22Tustison+N4ITK%22) papers and other relevant references in [Pubmed](http://www.ncbi.nlm.nih.gov/pubmed/?term=((%22Tustison+N%22)+AND+%22Avants+B%22)
*   Visualization: e.g. [a gource of ANTs development](http://vimeo.com/66781467)
*   [DiReCT](http://www.ncbi.nlm.nih.gov/pubmed/?term=((%22DIRECT%22)+AND+%22Avants%22+AND+DAS) cortical thickness [papers](http://www.ncbi.nlm.nih.gov/pubmed/?term=((%22Cortical+Thickness%22)+AND+%22Avants%22)
*   A [folder](https://sourceforge.net/projects/advants/files/Documentation/)
            of relevant docs: [segmentation](http://sourceforge.net/projects/advants/files/Documentation/atropos.pdf/download),
            [registration](http://sourceforge.net/projects/advants/files/Documentation/antstheory.pdf/download),
            [usage(old)](http://sourceforge.net/projects/advants/files/Documentation/ants.pdf/download),
            [for
            clinical apps](http://sourceforge.net/projects/advants/files/Documentation/ANTSMethodologySummary.docx/download)
*   ANTs redesigned for generality, automation,
                multi-core computation with ITKv4
*   Dev'd ITKv4 with Kitware, GE, Natl. Lib of
                Medicine & Academia
            </dd>

	    <dd class="clear"></dd>

            <dt>Acclaim for Robustness & Other News</dt>
            <dd>

## ANTs has won several competitions<span>
	         Unbiased & international</span>

*   New _antsRegistration_ [bash example](https://github.com/stnava/ANTs/blob/master/Scripts/newAntsExample.sh)
*   New _Eigenanatomy_ for [.nlm.nih.gov/pubmed/23269595">multivariate
	         neuroimage analysis](http://www.ncbi</a) via [PCA](http://www.ncbi.nlm.nih.gov/pubmed/23286132) &
	         [CCA](http://www.ncbi.nlm.nih.gov/pubmed/20083207) </a>
*   ANTs finished in 1st rank in [Klein
	         2009 intl. brain mapping competition](http://www.ncbi.nlm.nih.gov/pubmed/19195496)
*   ANTs finished 1st overall in [
	         EMPIRE10 intl. lung mapping competition](http://www.ncbi.nlm.nih.gov/pubmed/21632295)
*   ANTs is the standard registration for
	         [MICCAI-2013](http://www.miccai2013.org/) segmentation competitions
*   Conducting ANTs-based R tutorial @ MICCAI-2013
*   ITK-focused Frontiers in
	         Neuroinformatics research topic [here](http://www.frontiersin.org/neuroinformatics/researchtopics/neuroinformatics_with_the_insi/1580)
          </dd>

            <dd class="clear"></dd>

            <dt>Learning about ANTs</dt>
            <dd>
	      <li>Presentations: e.g. [a
	      Prezi about ANTs](http://prezi.com/mwrmcm-h9-w4/ants/?kw=view-mwrmcm-h9-w4&rc=ref-40024395) (WIP)</li>
	      <li>Reproducible science as a teaching tool: e.g. [ compilable ANTs tutorial](https://github.com/stnava/ANTS_MultiModality) (WIP)</li>
	      <li>Other examples [ slideshow](http://brianavants.wordpress.com) </li>
              <li>Landmark-based mapping for e.g. hippocampus [
            discussed here](https://sourceforge.net/p/advants/discussion/840261/thread/1cb7b165/?limit=50)</li>
	    <li>Brief ANTs segmentation [ video](http://vimeo.com/67814201)</li>
	    </dd>

            <dd class="clear"></dd>

            <dt>References</dt>
            <dd><li>[Google Scholar](http://scholar.google.com/scholar?q=Advanced+Normalization+Tools+%22ANTs%22+-ant&hl=en&as_sdt=1%2C39&as_ylo=2008&as_yhi=)</li>            
            <li>[Pubmed](http://www.ncbi.nlm.nih.gov/pubmed?term=(%22Avants%20B%22)%20OR%20%22Tustison%20N%22)</li>            
            <li>ANTs was supported by: R01-EB006266-01</li>
            <dd class="clear"></dd>
        </dl>

        <div class="clear"></div>

    </div>

</body>

</html>

