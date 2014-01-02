FILE=BP

all: $(FILE).pdf

$(FILE).pdf: *.tex Makefile
	make pdf
	bibtex $(FILE)
	@echo "Rerun to get" > $(FILE).log
	@while grep "Rerun to get" $(FILE).log >/dev/null 2>&1 ; do \
	make pdf; \
	done

pdf:
	pdflatex $(FILE).tex

clean:
	rm -f *~ *.aux *.log *.out *.toc *.bbl *.backup *.blg

distclean:
	make clean
	rm -f $(FILE).pdf
