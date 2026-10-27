# Run through tasks/01_ready.R, from the tutorial folder you downloaded or
# cloned (open qtl2mr-tutorial.Rproj first), in a fresh R session.
# Installs the bundled qtl2mr 0.6.0 into this project's private library/ folder.
cat("TASK 01 | Get ready\n")
archive <- "package/qtl2mr_0.6.0.tar.gz"
if (!file.exists(archive)) stop("Open qtl2mr-tutorial.Rproj first, then try again.")
if ("qtl2mr" %in% loadedNamespaces()) stop("Restart R, then run Task 01 again.")
# Check that the bundled archive is the unchanged qtl2mr 0.6.0 release file.
stopifnot(unname(tools::md5sum(archive)) == "163016881aaef1795df32541af460b99")
dir.create("library", showWarnings = FALSE)
workshop_lib <- normalizePath("library", winslash = "/")
.libPaths(c(workshop_lib, .libPaths()))
if (!requireNamespace("fs", quietly = TRUE)) {
  cat("Installing the fs dependency. This may need an internet connection.\n")
  install.packages("fs", lib = workshop_lib, repos = "https://cloud.r-project.org",
                   quiet = TRUE)
}
if (!requireNamespace("fs", quietly = TRUE)) stop("The fs installation needs attention.")
cat("Installing the supplied qtl2mr 0.6.0 archive...\n")
install.packages(archive, lib = workshop_lib, repos = NULL, type = "source",
                 quiet = TRUE)
library(qtl2mr, lib.loc = workshop_lib)
stopifnot(as.character(packageVersion("qtl2mr")) == "0.6.0",
          normalizePath(find.package("qtl2mr"), winslash = "/") ==
            normalizePath(file.path(workshop_lib, "qtl2mr"), winslash = "/"))
cat("qtl2mr version:", as.character(packageVersion("qtl2mr")), "\n")
