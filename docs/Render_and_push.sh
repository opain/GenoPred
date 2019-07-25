#!/bin/bash

#$ -cwd

module add general/R/3.5.0
module add general/pandoc/2.5

Rscript -e "rmarkdown::render_site()"
git add -u
git add .
git commit -m 'Update'
git push
