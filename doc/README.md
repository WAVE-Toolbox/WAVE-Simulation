# Documentation

- Complete documentation: `make all`

## Theory

- Make target `make install_theory`
- Select language by `LANGUAGE` to `DE` or `EN`
 - eg. `make install_theory LANGUAGE=EN` for english
- Output: `theory/WAVE_theory_${LANGUAGE}.pdf`
- Requirement: *pdflatex* and *bibtex*

## Documentation for source code

- Make target `make install_doxygen`
- Output: `doxygen/html/index.html`
- Requirement: *doxygen*
