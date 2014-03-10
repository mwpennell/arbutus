all: install

test: install
	make -C inst/tests test

# I dislike devtools' use of the Collate field (which causes problems
# for rapidly adding new code) so I'm disabling it this way:
DEVTOOLS_DOCUMENT=devtools::document(roclets=c('namespace', 'rd'))
# But Matt is using the Collate field so far, so let's leave it in:
DEVTOOLS_DOCUMENT=devtools::document()
document:
	@mkdir -p man
	Rscript -e "library(methods); ${DEVTOOLS_DOCUMENT}"

install:
	R CMD INSTALL --no-test-load .

build:
	R CMD build .

check: build
	R CMD check --no-manual `ls -1tr arbutus*gz | tail -n1`
	@rm -f `ls -1tr arbutus*gz | tail -n1`
	@rm -rf arbutus.Rcheck

# No real targets!
.PHONY: all test document install
