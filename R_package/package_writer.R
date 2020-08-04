# Let's build the primePCA package!

library(devtools)
library(roxygen2)
setwd("/Users/ziweizhu/Dropbox/BambooLab_Templates/R_package")
devtools::create("primePCA")
setwd("/Users/ziweizhu/Dropbox/BambooLab_Templates/R_package/primePCA") # In case the path is not switched to the package
devtools::document()
setwd("..")
devtools::install("primePCA")
