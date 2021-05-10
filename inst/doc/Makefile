## qs = 1 2 3 4 6 7 8 9 10 11 12 13 14 22 23 25 28

all: run

run:
	R --vanilla < Rmd-run.R

copy:
	cp q*.R /run/user/1000/gvfs/ftp\:host\=ftp.biostat3.net/download/R/
	cp index.html /run/user/1000/gvfs/ftp\:host\=ftp.biostat3.net/download/R/
	# cp installation.html /run/user/1000/gvfs/ftp\:host\=ftp.biostat3.net/download/R/
	cp labs/q*.html /run/user/1000/gvfs/ftp\:host\=ftp.biostat3.net/download/R/labs/
	cp solutions/q*.html /run/user/1000/gvfs/ftp\:host\=ftp.biostat3.net/download/R/solutions/

replace_once:
	cd labs && for file in 1 2 3 4 6 7 8 9 10 11 12 13 14 22 23 25 28; do cp q$$file.Rmd q$$file.Rmd~; sed -i '/output/,/---/coutput:\n  prettydoc::html_pretty:\n    theme: cayman\n    highlight: github\n---' q$$file.Rmd; done
	cd solutions && for file in 1 2 3 4 6 7 8 9 10 11 12 13 14 22 23 25 28; do cp q$$file.Rmd q$$file.Rmd~; sed -i '/output/,/---/coutput:\n  prettydoc::html_pretty:\n    theme: cayman\n    highlight: github\n---' q$$file.Rmd; done

