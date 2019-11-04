## qs = 1 2 3 4 6 7 8 9 10 11 12 13 22 23 25 28

all: copy_files

copy_files:
	sftp biostat3@biostat3.net lab/html

run:
	R --vanilla < Rmd-run.R

replace_once:
	cd labs && for file in 1 2 3 4 6 7 8 9 10 11 12 13 22 23 25 28; do cp q$$file.Rmd q$$file.Rmd~; sed -i '/output/,/---/coutput:\n  prettydoc::html_pretty:\n    theme: cayman\n    highlight: github\n---' q$$file.Rmd; done
	cd solutions && for file in 1 2 3 4 6 7 8 9 10 11 12 13 22 23 25 28; do cp q$$file.Rmd q$$file.Rmd~; sed -i '/output/,/---/coutput:\n  prettydoc::html_pretty:\n    theme: cayman\n    highlight: github\n---' q$$file.Rmd; done

