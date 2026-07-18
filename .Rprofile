reticulate::use_python("/usr/local/bin/python3", required = TRUE)

# Register knitr's figure-crop hook so `crop: true` (in _quarto.yml) actually
# trims whitespace around generated figures. hook_pdfcrop uses pdfcrop for PDF
# figures and ImageMagick `convert -trim` for PNG (matplotlib already crops via
# savefig.bbox=tight; this mainly tightens the R base plots' reserved title margin).
if (requireNamespace("knitr", quietly = TRUE)) {
  knitr::knit_hooks$set(crop = knitr::hook_pdfcrop)
}
