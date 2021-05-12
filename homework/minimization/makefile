CFLAGS    +=  $(shell gsl-config --cflags)
LDLIBS    +=  $(shell gsl-config --libs)
CC         = gcc
DISPLAY    = cat
LDLIBS    += -lm
OUTPUT     = higgsFit.txt
EXECS      = main
OBJECTS 	 = main.o minimization.o input2Array.o
HDRS    	 = minimization.h input2Array.h
FIGS    	 = higgsFit.png
TXT        = output.txt

.PHONEY: default
default: print $(FIGS)

print: $(TXT)
	$(DISPLAY) $<

higgsFit.png: higgsFit.gpi higgsFit.txt
	pyxplot $<

$(TXT): $(EXECS)
	./$< $(OUTPUT) 1> $@

$(EXECS): $(OBJECTS) $(HDRS)

clean:
	$(RM) $(EXECS) $(OBJECTS) $(TXT) $(OUTPUT) $(FIGS)
