## qs = 1 2 3 4 6 7 8 9 10 11 12 13 14 22 23 25 28

qs = 10 12 14 22 25 2 4 7 9 11 13 1 23 28 3 6 8
all: run

run:
	R --vanilla < Rmd-run.R

EXERCISE = 7
EXPR = 'library(rmarkdown); library(prettydoc); rmarkdown::render("q$(EXERCISE).Rmd", html_pretty())'

file:
	cd labs && Rscript -e $(EXPR)
	cd solutions && Rscript -e $(EXPR)

copy:
	cp q*.R /run/user/1000/gvfs/ftp\:host\=ftp.biostat3.net/download/R/
	cp index.html /run/user/1000/gvfs/ftp\:host\=ftp.biostat3.net/download/R/readme.html
	# cp installation.html /run/user/1000/gvfs/ftp\:host\=ftp.biostat3.net/download/R/
	cp labs/q*.html /run/user/1000/gvfs/ftp\:host\=ftp.biostat3.net/download/R/labs/
	cp solutions/q*.html /run/user/1000/gvfs/ftp\:host\=ftp.biostat3.net/download/R/solutions/

render:
	cd solutions && for ex in 10 12 14 22 25 2 4 7 9 11 13 1 23 28 3 6 8; do cp q$$ex.Rmd temp.Rmd && Rscript -e 'library(rmarkdown); library(prettydoc); rmarkdown::render("temp.Rmd", html_pretty())' && cp temp.html q$$ex.html ; done

replace_once:
	cd labs && for file in $qs; do cp q$$file.Rmd q$$file.Rmd~; sed -i '/output/,/---/coutput:\n  prettydoc::html_pretty:\n    theme: cayman\n    highlight: github\n---' q$$file.Rmd; done
	cd solutions && for file in $qs; do cp q$$file.Rmd q$$file.Rmd~; sed -i '/output/,/---/coutput:\n  prettydoc::html_pretty:\n    theme: cayman\n    highlight: github\n---' q$$file.Rmd; done

