# Let's build the primePCA package!

library(devtools)
library(roxygen2)
setwd("/Users/ziweizhu/Dropbox/BambooLab_Templates/R_package")
devtools::create("primePCA.test")
setwd("/Users/ziweizhu/Dropbox/BambooLab_Templates/R_package/primePCA.test") # In case the path is not switched to the package
# Please modify the DESCRIPTION file if your new package depends on other
# packages. Refer to the current DESCRIPTION file for details.  
devtools::document()
setwd("..")
devtools::install("primePCA_test")
