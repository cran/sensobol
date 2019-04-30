# sensobol 0.2.0

* New option in the `sobol_matrices` function: the option `cluster` allows to create Sobol' matrices for clusters of parameters.

* New test functions added: 
  - [Bratley & Fox 1988](https://dl.acm.org/citation.cfm?id=214372&dl=ACM&coll=DL)
  - [Bratley et al. 1992](https://dl.acm.org/citation.cfm?id=146385)
  - [Oakley & O'Hagan 2004](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.470.6932&rep=rep1&type=pdf)
  
* Added references to all test functions.

# sensobol 0.1.1

* The vignette rendered wrongly in the previous version;
now the issue is corrected.

* Corrected the following note, found in the CRAN check: 
  - Namespace in Imports field not imported from: ‘parallel’.
 
* Some functions that were exported to the R package manual 
are now internal.

# sensobol 0.1.0

* Added a NEWS.md file to track changes to the package.
