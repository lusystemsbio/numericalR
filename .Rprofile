# Point reticulate at a Python that actually has the book's packages (requirements.txt).
#
# Two clean-machine failures this guards against:
#  1. The path must NOT be hardcoded. This file used to open with
#     `use_python("/usr/local/bin/python3", required = TRUE)`, which is where the
#     python.org installer puts python3 but NOT where Homebrew does (/opt/homebrew/bin on
#     Apple Silicon). When that path is absent, `required = TRUE` throws, and because it
#     was the FIRST line it aborted this whole file, so the knitr crop/rmar hooks below
#     never registered either.
#  2. With no Python configured, reticulate falls back to a "managed" ephemeral venv: it
#     downloads uv plus its own CPython, installs only the packages it can infer, and the
#     render then dies mid-chapter on e.g. "ModuleNotFoundError: No module named 'pandas'"
#     while the real interpreter sat there fully provisioned. RETICULATE_USE_MANAGED_VENV
#     turns that off, so a bad lookup fails loudly instead of silently substituting.
Sys.setenv(RETICULATE_USE_MANAGED_VENV = "no")

local({
  # A set RETICULATE_PYTHON overrides use_python() (reticulate warns that the call "will
  # be ignored"), so honor an explicit override and let reticulate report it if it is wrong.
  if (nzchar(Sys.getenv("RETICULATE_PYTHON"))) return(invisible(NULL))

  needed <- c("numpy", "scipy", "matplotlib", "pandas")
  cands <- c("/usr/local/bin/python3",                  # python.org installer
             "/opt/homebrew/bin/python3",               # Homebrew (Apple Silicon)
             "/usr/local/opt/python3/bin/python3",      # Homebrew (Intel)
             Sys.which("python3"),
             rev(sort(Sys.glob("/Library/Frameworks/Python.framework/Versions/*/bin/python3"))))
  cands <- unique(cands[nzchar(cands) & file.exists(cands)])

  # Accept the first candidate that can actually import the book's packages, so an
  # interpreter that merely exists but is unprovisioned is skipped rather than chosen.
  probe <- sprintf("import %s", paste(needed, collapse = ", "))
  for (py in cands) {
    if (identical(0L, suppressWarnings(system2(py, c("-c", shQuote(probe)),
                                               stdout = FALSE, stderr = FALSE)))) {
      reticulate::use_python(py, required = TRUE)
      return(invisible(NULL))
    }
  }
  message("[.Rprofile] No python3 found with the book's packages (",
          paste(needed, collapse = ", "), ").\n",
          "  Checked: ", if (length(cands)) paste(cands, collapse = ", ") else "no python3 at all",
          "\n  Fix: python3 -m pip install -r requirements.txt",
          "\n  or set RETICULATE_PYTHON to the interpreter that has them.")
})

# Register knitr's figure-crop hook so `crop: true` (in _quarto.yml) actually
# trims whitespace around generated figures. hook_pdfcrop uses pdfcrop for PDF
# figures and ImageMagick `convert -trim` for PNG (matplotlib already crops via
# savefig.bbox=tight; this mainly tightens the R base plots' reserved title margin).
if (requireNamespace("knitr", quietly = TRUE)) {
  knitr::knit_hooks$set(crop = knitr::hook_pdfcrop)

  # R base plots reserve a wide top margin for a title (usually absent here). The
  # crop hook trims it, which leaves the figure WIDE (~7:4.2) rather than the 7:5
  # matplotlib produces, so R and Python plots came out different sizes. Shrink the
  # top/right margin (enabled per chunk via `rmar: true` in _quarto.yml opts_chunk)
  # so a cropped R landscape plot keeps ~7:5. Chunks that set their own par(mar=...)
  # or par(pty="s") still override/extend this default.
  knitr::knit_hooks$set(rmar = function(before, options) {
    if (before && isTRUE(options$rmar) && grepl("^r", tolower(options$engine))) {
      par(mar = c(4.6, 4.6, 0.5, 1.0))   # top=0.5 -> cropped landscape ~1.40, matching matplotlib
    }
  })
}
