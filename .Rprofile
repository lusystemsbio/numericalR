reticulate::use_python("/usr/local/bin/python3", required = TRUE)

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
      par(mar = c(4.6, 4.6, 1.0, 1.0))
    }
  })
}
