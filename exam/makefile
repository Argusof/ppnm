CFLAGS    +=  $(shell gsl-config --cflags)
LDLIBS    +=  $(shell gsl-config --libs)
CC         = gcc
DISPLAY    = cat
LDLIBS    += -lm
OUTPUT     = outputFile.txt
EXECS      = main
OBJECTS 	 = main.o jacobi.o utilities.o
HDRS    	 = jacobi.h utilities.h
FIGS    	 =
TXT        = output.txt

.PHONEY: default
default: print #$(FIGS)

print: $(TXT)
	#$(DISPLAY) $<

$(TXT): $(EXECS)
	./$< $(OUTPUT) #1> $@

$(EXECS): $(OBJECTS) $(HDRS)

clean:
	$(RM) $(EXECS) $(OBJECTS) $(TXT) $(OUTPUT) $(FIGS)
