require(knitr)
require(markdown)
rmarkdown::render('matrixpls-intro.Rmd', clean = FALSE, run_pandoc = FALSE)

system("/Applications/RStudio.app/Contents/MacOS/pandoc/pandoc +RTS -K512m -RTS matrixpls-intro.utf8.md --to latex --from markdown+autolink_bare_uris+ascii_identifiers+tex_math_single_backslash --output matrixpls-intro.tex --template /Library/Frameworks/R.framework/Versions/3.3/Resources/library/rticles/rmarkdown/templates/jss_article/resources/template.tex --highlight-style tango --latex-engine /opt/local/bin/pdflatex --bibliography matrixpls-intro.bib --filter /Applications/RStudio.app/Contents/MacOS/pandoc/pandoc-citeproc")
