check:
	R CMD build .
	R CMD check `ls *.tar.gz`
