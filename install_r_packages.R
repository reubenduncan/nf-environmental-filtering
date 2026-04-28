# install_r_packages.R
# Installs CRAN packages not available on conda-forge into the active conda
# environment's R library.
#
# Run once on the login node after `conda env create`:
#   conda run -n environmental-filtering Rscript \
#       projects/environmental-filtering/install_r_packages.R

pkgs <- c("NST", "ecodist", "metacom")

conda_prefix <- Sys.getenv("CONDA_PREFIX")
if (nchar(conda_prefix) == 0) {
  stop("CONDA_PREFIX is not set. Activate the conda environment first, or use:\n",
       "  conda run -n environmental-filtering Rscript install_r_packages.R")
}

lib <- file.path(conda_prefix, "lib", "R", "library")
if (!dir.exists(lib)) {
  stop("R library directory not found: ", lib,
       "\nIs R installed in the conda environment?")
}

message("Installing into: ", lib)

for (pkg in pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    message(pkg, ": already installed, skipping.")
    next
  }
  message("Installing ", pkg, " ...")
  install.packages(pkg, lib = lib, repos = "https://cloud.r-project.org")
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Installation of '", pkg, "' failed. Check warnings above.")
  }
  message(pkg, ": installed successfully.")
}

message("All required packages are installed.")
