CFLAGS    +=  $(shell gsl-config --cflags)
LDLIBS    +=  $(shell gsl-config --libs)
CC         = gcc
DISPLAY    = cat
LDLIBS    += -lm
OUTPUT     = outputFile.txt #myTime.txt gslTime.txt
EXECS      = main
OBJECTS 	 = main.o integration.o #utilities.o diffClock.o test_runtime.o
HDRS    	 = integration.h #utilities.h diffClock.h
FIGS    	 = #energies.png #myTime.png gslTime.png
TXT        = output.txt

.PHONEY: default
default: print #$(FIGS)

print: $(TXT)
	$(DISPLAY) $<

#energies.png: energies.gpi outputFile.txt
#	pyxplot $<

#myTime.png: myTime.gpi myTime.txt
#	pyxplot $<

#energies.png: gslTime.gpi gslTime.txt
#	pyxplot $<

$(TXT): $(EXECS)
	./$< $(OUTPUT) 1> $@

$(EXECS): $(OBJECTS) $(HDRS)

clean:
	$(RM) $(EXECS) $(OBJECTS) $(TXT) $(OUTPUT) $(FIGS)
