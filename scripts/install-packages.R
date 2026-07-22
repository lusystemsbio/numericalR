# R packages required to render this book. Run once after cloning:
#
#   Rscript scripts/install-packages.R
#
# Installs the latest CRAN version of each; versions confirmed working as of
# this script's last update are noted alongside each package. If a chapter
# starts failing to render and a package version is suspected, cross-check
# against these first.

pkgs <- c(
  "reticulate",      # 1.46.0 -- R/Python bridge (see .Rprofile, _quarto.yml execute-dir)
  "languageserver",  # 0.3.18 -- VS Code R language support
  "deSolve",         # 1.42   -- ODE/DDE integration (Part 2 onward)
  "microbenchmark",  # 1.5.0  -- timing comparisons
  "ggplot2",         # 4.0.3
  "rootSolve",       # 1.8.2.4 -- root finding (2E/2F)
  "magick",          # 2.9.1  -- enables the crop: true figure-trim hook
  "patchwork",       # 1.3.2  -- combining ggplots side-by-side (9A)
  "Rtsne",           # 0.17   -- t-SNE (10A/10B)
  "uwot",            # 0.2.4  -- UMAP (10A/10B)
  "mclust",          # 6.1.3  -- Gaussian mixture models (10B)
  "igraph",          # 2.3.3  -- graph drawing, betweenness, community detection (10C)
  "princurve",       # 2.1.6  -- principal curve (10A.4)
  "bio3d"            # 2.4.5  -- PDB structure parsing (1D.1)
)

installed <- rownames(installed.packages())
missing <- setdiff(pkgs, installed)
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  cat("All required packages already installed.\n")
}
