# Renders the report in both formats:
#   phage_datamining.html  (HTML)
#   phage_datamining.md    (readable directly on GitHub, figures in figures/)
rmarkdown::render("phage_datamining.Rmd", output_format = "all")

