
wd <- getwd()
req_lib_dir <- paste0(wd,"/req_lib")
packageVersion("rlang")


###################### Creating folders:
subfolder_new        <- paste0("req_lib/")

if (!dir.exists(subfolder_new)) {
       dir.create(subfolder_new)
}




if (!require(mvtnorm, lib = req_lib_dir)) {
    .libPaths(req_lib_dir)
    install.packages("mvtnorm", lib = req_lib_dir,
                     repos = "https://archive.linux.duke.edu/cran/")
    library(mvtnorm, lib = req_lib_dir)
}

if (!require(SimDesign, lib = req_lib_dir)) {
    .libPaths(req_lib_dir)
    install.packages("SimDesign", lib = req_lib_dir,
                     repos = "https://archive.linux.duke.edu/cran/")
    library(SimDesign, lib = req_lib_dir)
}

if(!require(MASS, lib = req_lib_dir)) {
    .libPaths(req_lib_dir)
    install.packages("MASS", lib = req_lib_dir,
                     repos = "https://archive.linux.duke.edu/cran/")
    library(MASS, lib = req_lib_dir)
}

if (!require(equalCovs, lib = req_lib_dir)) {
    .libPaths(req_lib_dir)
    install.packages(
        'equalCovs_1.0.tar.gz',
        repos=NULL, type='source',
        lib = req_lib_dir)
    library(equalCovs, lib = req_lib_dir)
}