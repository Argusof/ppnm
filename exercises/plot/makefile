CFLAGS = -O -Wall -std=gnu11 $(shell gsl-config --cflags)
LDLIBS = -lm $(shell gsl-config --libs)
FIGS = pyx_erf.png pyx_gamma.png
EXCES = main
OBJECTS = main.o my_erf.o my_gamma.o
INPUT = input_erf.txt input_gamma.txt
OUTPUT = output_erf.txt output_gamma.txt
TMP = tmp_erf.gpi tmp_gamma.gpi
HDRS = my_erf.h my_gamma.h

.PHONEY:default
default: $(FIGS)

pyx_erf.png: output_erf.txt input_erf.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key bottom right;\
	set tics in;\
	set xlabel "x";\
	set ylabel "erf(x)";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line lw 2 title "erf from math.h"\
	,"$<" using 1:3 with line lw 2 title "erf from gsl"\
	,"$<" using 1:4 with line lw 2 title "my erf"\
	,"input_erf.txt" using 1:2 with points pointtype 2 title "tabulated data"\
	' | pyxplot

pyx_gamma.png: output_gamma.txt input_gamma.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key bottom right;\
	set tics in;\
	set xlabel "x";\
	set ylabel "gamma(x)";\
	set title "Gamma function";\
	plot \
	 "$<" using 1:2 with line lw 2 title "gamma from math.h"\
	,"$<" using 1:3 with line lw 2 title "gamma from gsl"\
	,"$<" using 1:4 with line lw 2 title "my gamma"\
	,"input_gamma.txt" using 1:2 with points pointtype 2 title "tabulated data"\
	' | pyxplot

$(TMP): pyx_erf.png pyx_gamma.png
		touch $@

$(OUTPUT): $(EXCES)
		./$< output_erf.txt
		./$< output_gamma.txt

$(EXCES): $(OBJECTS)

$(OBJECTS): $(HDRS)

.PHONEY:clean
clean:
		$(RM) $(EXCES) $(OBJECTS) $(FIGS) $(OUTPUT) $(TMP)
