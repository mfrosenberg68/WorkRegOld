all: mplib.PS

mplib.dvi: mplib.tex mplib.bbl
	latex mplib

clean:
	rm -f mplib.dvi mplib.log

.SUFFIXES:
.SUFFIXES: .dvi .PS $(SUFFIXES)

.dvi.PS:
	dvips -o $@ $*
