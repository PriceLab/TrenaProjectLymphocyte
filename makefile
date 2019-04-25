all:  docs install

docs:
	R -e "devtools::document()"
build:
	(cd ..; R CMD build --no-build-vignettes TrenaProjectLymphocyte)

install:
	(cd ..; R CMD INSTALL TrenaProjectLymphocyte)

check:
	(cd ..; R CMD check `ls -t TrenaProjectLymphocyte) | head -1`)

test:
	for x in inst/unitTests/test_*.R; do echo $$x; R -f $$x; done

tv:
	R -f inst/demos/trenaViz/tv.R 
