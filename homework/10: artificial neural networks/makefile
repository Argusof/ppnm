CFLAGS    +=  $(shell gsl-config --cflags)
LDLIBS    +=  $(shell gsl-config --libs)
CC         = gcc
DISPLAY    = cat
LDLIBS    += -lm
INPUT      = cos_tabulated_20.txt data.txt
EXECS      = main
OBJECTS 	 = main.o minimization.o neuralNetwork.o utilities.o input2Array.o
HDRS    	 = minimization.h neuralNetwork.h utilities.h input2Array.h
FIGS    	 = plot.png
TXT        = output.txt

.PHONEY: default
default: print $(FIGS)

print: $(TXT)
	$(DISPLAY) $<

plot.png: plot.gpi data.txt
	pyxplot $<

$(TXT): $(EXECS)
	./$< $(INPUT) 1> $@

$(EXECS): $(OBJECTS) $(HDRS)

clean:
	$(RM) $(EXECS) $(OBJECTS) $(TXT) $(OUTPUT) $(FIGS)
