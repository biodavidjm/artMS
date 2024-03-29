# artMS

___Analytical R Tools for Mass Spectrometry___

---

[![Build Status](https://travis-ci.com/biodavidjm/artMS.svg?branch=master)](https://travis-ci.com/biodavidjm/artMS)
[![codecov](https://codecov.io/github/biodavidjm/artMS/branch/master/graphs/badge.svg)](https://codecov.io/github/biodavidjm/artMS) 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8093247.svg)](https://doi.org/10.5281/zenodo.8093247)

## Overview

`artMS` ([http://artms.org/](http://artms.org/)) is an R package that provides a set of tools for the analysis and integration of large-scale proteomics (mass-spectrometry-based) datasets obtained using the popular proteomics software package 
[MaxQuant](http://www.biochem.mpg.de/5111795/maxquant). The [functions available in artMS](https://biodavidjm.github.io/artMS/reference/index.html) can be grouped into 4 major categories:

- Multiple quality control (QC).
- Relative quantification using [MSstats](http://msstats.org/).
- Downstream analysis and integration of quantifications (enrichment, clustering, PCA, summary plots, etc)
- Generation of input files for other tools, including [SAINTq and SAINTexpress](http://saint-apms.sourceforge.net/Main.html), [Photon](https://github.com/jdrudolph/photon), and [Phosfate](http://phosfate.com/)


`artMS` performs the different analyses taking as input the following files:

- `evidence.txt` file: The output of the quantitative proteomics software 
package `MaxQuant`. 
- `keys.txt` (tab-delimited) txt file generated by the user describing the experimental designed (check below to learn how to create it).
- `contrast.txt` (tab-delimited) txt file generated by the user with the comparisons between conditions to be quantified (check below to learn how to create it).
- `config.yaml`: a configuration file which enables the customization of a number of parameters for the quantification (and other operations, including QC analyses, charts and annotations). A configuration file template can be generated by running `artmsWriteConfigYamlFile()`



## How to install

### Bioconductor

`artMS version >= 1.10.1` had many changes to adjust for changes in MSstats. This version requires:

- Install `R version >= 4.1.0` (check the R version running on your system by executing the function `getRversion()`)
- Bioconductor: `BiocManager::install("BiocVersion")`
- artMS: `BiocManager::install("artMS")`
- If you are planning to use the `artmsAnalysisQuantifications()` to perform a comprehensive downstream analysis of the quantitative results, then install the following packages:

```
# From bioconductor:
BiocManager::install(c("ComplexHeatmap", "org.Mm.eg.db"))

# From CRAN:
install.packages(c("factoextra", "FactoMineR", "gProfileR", "PerformanceAnalytics"))
```

Extra: Why Bioconductor? [Here you can find a nice summary of many good reasons](https://bioinformatics.stackexchange.com/questions/639/why-bioconductor)).

### Development version from Github (unstable)

Assuming that you have an `R (>= 4.1)` version running on your system, 
follow these steps:

```
install.packages("devtools")
library(devtools)
install_github("biodavidjm/artMS")
```

Once installed, the package can be loaded and attached to your current 
workspace as follows:

```{r, eval=TRUE}
library(artMS)
```

Once installed, we suggest you to do a quick test by running the quality control functions using the "evidence" (`artms_data_ph_evidence`) and "keys" (`artms_data_ph_keys`) files included in `artMS` as test datasets.

```
# First go to a local working directory: several pdfs will be generated
# setwd("/path/to/your/working/directory/")

# And run:
artmsQualityControlEvidenceBasic(evidence_file = artms_data_ph_evidence,
                                  keys_file = artms_data_ph_keys, 
                                  prot_exp =  "PH")
```

(To learn more about these testing datasets, check the documentation by running `?artms_data_ph_keys` or `?artms_data_ph_evidence` on the R console)


Once the QC is done, go to the folder `"/path/to/your/working/directory/"` and check out all the generated QC (pdf) files available in the `qc_basic` folder

## How to Contribute to artMS

`artMS` is an open source project, therefore you are more than welcome to contribute and make the analysis of Mass Spectrometry data easier and better using this fantastic language and environment for statistical computing and graphics (i.e. `R`).

There are multiple options:

- [Submit issues to this repo](https://github.com/biodavidjm/artMS/issues) reporting problems, bugs, or suggesting new features.
- Fork and make pull requests. To find out more about this option, 
some very useful guides for beginners can be found <a href="https://akrabat.com/the-beginners-guide-to-contributing-to-a-github-project/" target="blank">here</a>
and <a href="https://github.com/Bioconductor/Contributions/blob/master/CONTRIBUTING.md" target="blank">there</a> (or even <a href="http://lmgtfy.com/?q=how+to+contribute+to+a+github+project" target="blank">beyond</a>). When submitting a Pull Request, don't forget to select @biodavidjm as the reviewer

__Tips__: Do you need to remember the basics of markdown? [Check out this fantastic link](https://commonmark.org/help/tutorial/index.html).


## artMS Help available online

- The vignette can also be accessed at [http://artms.org](http://artms.org)
- Errors or warnings? Please, 
<a href="https://github.com/biodavidjm/artMS/issues" target="_blank">submit them as a new issue</a>
at the official Github repository
- Any other inquiries: <artms.help@gmail.com>





