#!/bin/bash
doc=expose
latex $doc.tex > poubelle
dvips -q $doc.dvi > poubelle
ps2pdf $doc.ps > poubelle
#mv $doc.pdf expose.pdf
rm poubelle
