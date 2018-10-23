check:
	R-devel CMD build .
	R-devel CMD check --as-cran biostat3_`grep Version DESCRIPTION | cut -b 10-15`.tar.gz
