LANG=DE

ifeq ($(LANG),DE)
	TARGET=WAVE_theory_DE
else
	TARGET=WAVE_theory_EN
endif

.PHONY: install
install: pdf

.PHONY: all
all: pdf

pdf:
	pdflatex $(TARGET)
	bibtex $(TARGET)
	pdflatex $(TARGET)
	pdflatex $(TARGET)
	make clean

.PHONY: pdf
clean:
	rm -rf *.aux *.lof *.log *.lot *.toc *.bbl *.blg *.out
