# Makefile for Sphinx documentation
#

SHELL := /bin/bash

# You can set these variables from the command line.
SPHINXOPTS    = -j8 -v
SPHINXBUILD   ?= LANG=C sphinx-build
PAPER         =

BUILD_DIR     ?= build
FILES         =

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILD_DIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) source

.PHONY: help clean html roadmap

#------------------------------------------------------------------------------
all: html
help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html      to make standalone HTML files"
	@echo "  latex     to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  roadmap   to regenerate source/roadmap.rst from changelog.rst"

clean:
	-rm -rf $(BUILD_DIR)/* source/generated source/api/generated

roadmap:
	python3 generate_roadmap.py

html: roadmap
	mkdir -p $(BUILD_DIR)/html $(BUILD_DIR)/doctrees
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILD_DIR)/html $(FILES)
	@echo
	@echo "Build finished. The HTML pages are in build/html."

latex:
	mkdir -p build/latex build/doctrees
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) build/latex $(FILES)

html_check:
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) build/html $(FILES)

