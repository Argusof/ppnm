CFLAGS    +=  $(shell gsl-config --cflags)
LDLIBS    +=  $(shell gsl-config --libs)
CC         = gcc
DISPLAY    = cat
LDLIBS    += -lm
OUTPUT     = hydrogen.txt convData.txt
EXECS      = main
OBJECTS 	 = main.o rootfinding.o GSdecomp.o backsub.o rkstep.o
HDRS    	 = rootfinding.h GSdecomp.h backsub.h rkstep.h
FIGS    	 = hydrogen.png conv.png
TXT        = output.txt

.PHONEY: default
default: print $(FIGS)

print: $(TXT)
	$(DISPLAY) $<

hydrogen.png: hydrogen.gpi hydrogen.txt
	pyxplot $<

conv.png: conv.gpi convData.txt
	pyxplot $<

$(TXT): $(EXECS)
	./$< $(OUTPUT) 1> $@

$(EXECS): $(OBJECTS) $(HDRS)

clean:
	$(RM) $(EXECS) $(OBJECTS) $(TXT) $(OUTPUT) $(FIGS)
