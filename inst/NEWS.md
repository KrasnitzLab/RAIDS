CHANGES IN VERSION 1.5.1
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o Integration of Rsamtools

CHANGES IN VERSION 1.3.3
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o More comprehensive vignette

CHANGES IN VERSION 1.3.2
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o New functions inferAncestry(), inferAncestryGeneAware() and getRefSuperPop() to simplify ancestry inference

CHANGES IN VERSION 1.3.1
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o New function createAccuracyGraph() that creates a graphic representation of the accuracy for different values of PCA dimensions and K-neighbors through all tested ancestries.

CHANGES IN VERSION 0.99.15
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o Updating installation section in vignette.

CHANGES IN VERSION 0.99.14
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o Adding missing author David Tuveson.
    o Updating BiocViews terms.
    
CHANGES IN VERSION 0.99.13
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o Update in Reference GDS vignette.

CHANGES IN VERSION 0.99.12
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o Seven new loadable objects are available in the package.
    o The new readSNVVCF() function enable the use of VCF SNP files as input for the runExomeAncestry() and runRNAAncestry() functions.
    

CHANGES IN VERSION 0.99.11
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o Update main vignette.

CHANGES IN VERSION 0.99.10
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o Update main vignette.

CHANGES IN VERSION 0.99.9
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o Better documentation for the runRNAAncestry() function.


CHANGES IN VERSION 0.99.8
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o Some examples have been updated in the documentation.


CHANGES IN VERSION 0.99.7
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o The new runRNAAncestry() function executes most steps leading to the ancestry inference call on a specific RNA profile.
    o A vignette describing the content of the Reference GDS files has been created.
    o More parameter names have been changed to follow the camelCase style.
    

CHANGES IN VERSION 0.99.6
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o The vignette now referes to the generic formatted reference GDS rather than 1KG GDS file to showcase that the software is not dependant of the 1KG GDS file. Any refence dataset can be used as long as the dataset is formatted into a GDS file.
    

CHANGES IN VERSION 0.99.5
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o The function documentation has been improved.
    o New vignette has been created. The vignette covers the steps done by the runExomeAncestry() function.
    

CHANGES IN VERSION 0.99.4
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o More parameter names have been changed to follow the camelCase style.
    o The function documentation has been improved.
    o A wrapper function runExomeAncestry() is now available.
    

CHANGES IN VERSION 0.99.3
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o More parameter names have been changed to follow the camelCase style.
    
BUG FIXES

    o The warning related to the package man page has been removed. 
    

CHANGES IN VERSION 0.99.2
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o Most parameter names have been changed to follow the camelCase style.
    

CHANGES IN VERSION 0.99.1
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o The new runExomeAncestry() function encapsulates multiple ancestry inference steps in one command.
    
BUG FIXES

    o Ensure GDS file is closed before using stop() in the addPhase1KG2SampleGDSFromFile() function. 

CHANGES IN VERSION 0.99.0
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o The number of visible functions has been limited to simplify usage.

CHANGES IN VERSION 0.50.1
------------------------

SIGNIFICANT USER-VISIBLE CHANGES

    o NEWS file added to the package
