INFILE=main.tex
LATEX=pdflatex

# TODO : svg and png support

all:
	$(LATEX) $(INFILE)
	bibtex $(INFILE:.tex=)
	$(LATEX) $(INFILE)

clean:
	rm -f *.pdf *.out *aux *bbl *blg *log *toc *.ptb *.tod *.fls *.fdb_latexmk *.lof