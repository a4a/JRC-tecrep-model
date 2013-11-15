

library(colorout)
library(knitr)

purl("EWG-12-19-a4a-content.Rnw")

knit("EWG-12-19-a4a-content.Rnw")

system("pdflatex -synctex=1 -interaction=nonstopmode EWG-12-19-a4a-content.tex")



