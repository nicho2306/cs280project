paper.pdf: paper.tex *.tex *.bib *.sty $(wildcard figures/*)
	pdflatex $(patsubst %.tex,%,$<)
	bibtex $(patsubst %.tex,%,$<)
	pdflatex $(patsubst %.tex,%,$<)
	pdflatex $(patsubst %.tex,%,$<)

clean:
	rm -f *.aux *.bbl *.blg *.log *.pdf *.out
