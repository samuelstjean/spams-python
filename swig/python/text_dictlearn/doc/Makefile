OPT=-O article.hva macrocss.hva
#fancysection.hva
all:
#	pdflatex doc_dictlearn.tex
#	bibtex doc_dictlearn
	pdflatex doc_dictlearn.tex
	pdflatex doc_dictlearn.tex
hevea:
	hevea ${OPT} doc_dictlearn.tex
#	bibhva doc_dictlearn
	hevea ${OPT} doc_dictlearn.tex
	hevea ${OPT} doc_dictlearn.tex
	hacha -tocbis doc_dictlearn.html
	cp style.css doc_dictlearn.css
	mv *.html html/
	mv *.gif html/
	mv *.css html/
	mv html/style.css ./
clean:
	rm -f *.aux *.bbl *.haux *.hbbl *.blg *.toc *.htoc *.log

