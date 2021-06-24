CFLAGS    +=  $(shell gsl-config --cflags)
LDLIBS    +=  $(shell gsl-config --libs)
CC         = gcc
DISPLAY    = cat
LDLIBS    += -lm
INPUT      = data.txt
OUTPUT     = fit.plt
EXECS      = main
OBJECTS 	 = main.o GSdecomp.o backsub.o fillArrayFromInputData.o leastSquares.o
HDRS    	 = GSdecomp.h backsub.h leastSquares.h
FIGS    	 = fit.png
TXT        = output.txt

.PHONEY: default
default: print $(FIGS)

print: $(TXT)
	$(DISPLAY) $<

fit.png: fit.gpi $(INPUT) $(OUTPUT)
	pyxplot $<

$(TXT): $(EXECS) $(INPUT)
	./$< $(INPUT) $(OUTPUT) 1> $(TXT)

$(EXECS): $(OBJECTS) $(HDRS)

clean:
	$(RM) $(EXECS) $(OBJECTS) $(FIGS) $(TXT) $(OUTPUT)
