# rankrate 1.2.1

* Minor fixes to package vignettes to meet CRAN policy.

# rankrate 1.2.0

* Removed bugs in estimation of object quality vector, p.
* Removed "LP" estimation option, which is overly complex and has no apparent estimation advantages. We recommend using the "ASTAR" option instead, which provides exact MLE search and is generally faster.
* Allow for identification of multiple MLEs, if ties between objects. When this occurs, estimation functions produce a message.
* Added vignettes.

# rankrate 1.1.0

* First release.
