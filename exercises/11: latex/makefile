CFLAGS = -Wall -O1 -std=gnu11
LDLIBS = -lm

theExpFunc.pdf: theExpFunc.tex myExpFunc-fig-pyxplot.pdf
	pdflatex $<
	pdflatex $<

out_data.txt:main
	./$< > $@

main:main.o myExp.o

myExpFunc-fig-pyxplot.pdf: myExpFunc-fig-pyxplot.gpi out_data.txt
	pyxplot $<

clean:
	$(RM) theExpFunc.pdf *.log *.aux out* main *.o fig*.pdf fig*.tex

