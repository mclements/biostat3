check: build
	R-devel CMD check --as-cran biostat3_`awk '/Version:/ {print $$2}' DESCRIPTION`.tar.gz

build:
	R-devel CMD build .

