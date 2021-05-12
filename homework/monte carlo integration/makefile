CFLAGS    +=  $(shell gsl-config --cflags)
LDLIBS    +=  $(shell gsl-config --libs)
CC         = gcc
DISPLAY    = cat
LDLIBS    += -lm
OUTPUT     = errors.txt
EXECS      = main
OBJECTS 	 = main.o MCintegration.o
HDRS    	 = MCintegration.h
FIGS    	 = errorScaling.png
TXT        = errors.txt

.PHONEY: default
default: print $(FIGS)

print: $(TXT)
	$(DISPLAY) $<

errorScaling.png: errors.gpi errors.txt
	pyxplot $<

$(TXT): $(EXECS)
	./$< $(OUTPUT) #1> $@

$(EXECS): $(OBJECTS) $(HDRS)

clean:
	$(RM) $(EXECS) $(OBJECTS) $(TXT) $(OUTPUT) $(FIGS)
