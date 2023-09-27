<!-- badges: start -->
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![codecov](https://codecov.io/gh/krasnitzlab/RAIDS/branch/main/graph/badge.svg?token=LPFLOMUDVT)](https://codecov.io/gh/krasnitzlab/RAIDS)
[![R-CMD-check-bioc](https://github.com/krasnitzlab/RAIDS/actions/workflows/check-bioc.yaml/badge.svg)](https://github.com/krasnitzlab/RAIDS/actions/workflows/check-bioc.yaml)
<!-- badges: end -->


# Accurate genetic ancestry inference from cancer-derived molecular data with **RAIDS** #

The Robust Ancestry Inference using Data Synthesis (**_RAIDS_**) package 
enables accurate and robust inference of genetic 
ancestry from various types of molecular data, including whole-genome, 
whole-exome, targeted gene panels and RNA sequences, as described in our 
manuscript. Our tools retain high accuracy in presence of somatic 
alterations, such as those caused by cancer.

**This code and analysis pipeline were designed and developed for the following publication:**

>  Pascal Belleau, Astrid Deschênes, Nyasha Chambwe, David A. Tuveson, Alexander Krasnitz; Genetic Ancestry Inference from Cancer-Derived Molecular Data across Genomic and Transcriptomic Platforms. Cancer Res 1 January 2023; 83 (1): 49–58. https://doi.org/10.1158/0008-5472.CAN-22-0682


## Authors ##

[Pascal Belleau](http://ca.linkedin.com/in/pascalbelleau "Pascal Belleau"),
[Astrid Desch&ecirc;nes](https://www.linkedin.com/in/astriddeschenes "Astrid Desch&ecirc;nes"),
[David A. Tuveson](https://tuvesonlab.labsites.cshl.edu/) and
[Alexander Krasnitz](https://www.cshl.edu/research/faculty-staff/alexander-krasnitz/ "Alexander Krasnitz")


## Citing ##

If you use the **_RAIDS_** package for a publication, we would ask you to cite 
the following:

>  Pascal Belleau, Astrid Deschênes, Nyasha Chambwe, David A. Tuveson, Alexander Krasnitz; Genetic Ancestry Inference from Cancer-Derived Molecular Data across Genomic and Transcriptomic Platforms. Cancer Res 1 January 2023; 83 (1): 49–58. https://doi.org/10.1158/0008-5472.CAN-22-0682


## Installation ##

To install the latest version accessible, the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) 
package is required.

     ## Load required package
     library(devtools)

     ## Install the latest version of RAIDS
     devtools::install_github('KrasnitzLab/RAIDS')


## Documentation ##

[RAIDS Website](https://krasnitzlab.github.io/RAIDS/)

[RAIDS Get Started](https://krasnitzlab.github.io/RAIDS/articles/RAIDS.html)


## License ##

This package and the underlying **_RAIDS_** code are distributed under 
the Apache-2.0 license. You are free to use and redistribute this software. 

For more information on Apache-2.0 License see
[https://opensource.org/licenses/Apache-2.0](https://opensource.org/licenses/Apache-2.0)



## Maintainer

[Pascal Belleau](https://github.com/belleau/ "Pascal Belleau")


## Bugs/Feature requests ##

[Please contact us](https://github.com/KrasnitzLab/RAIDS/issues) for bug fixes or with feature requests. 

Thanks!
