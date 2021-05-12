CFLAGS    +=  $(shell gsl-config --cflags)
LDLIBS    +=  $(shell gsl-config --libs)
CC         = gcc
DISPLAY    = cat
LDLIBS    += -lm
OUTPUT     = harmOutputFile.txt SIRoutputFile.txt SIRoutputFile2.txt threeBodyOutputFile.txt
EXECS      = main
OBJECTS 	 = main.o rkstep.o
HDRS    	 = rkstep.h
FIGS    	 = harmonic.png SIR.png SIR2.png threeBody.png
TXT        = output.txt

.PHONEY: default
default: print $(FIGS)

print: $(TXT)
	$(DISPLAY) $<

#$(FIGS): harmonic.gpi SIR.gpi $(OUTPUT)
#	pyxplot $<

harmonic.png: harmonic.gpi harmOutputFile.txt
	pyxplot $<

SIR.png: SIR.gpi SIRoutputFile.txt
	pyxplot $<

SIR2.png: SIR2.gpi SIRoutputFile2.txt
	pyxplot $<

threeBody.png: threeBody.gpi threeBodyOutputFile.txt
	pyxplot $<

$(TXT): $(EXECS)
	./$< $(OUTPUT) 1> $(TXT)

$(EXECS): $(OBJECTS) $(HDRS)

clean:
	$(RM) $(EXECS) $(OBJECTS) $(TXT) $(OUTPUT) $(FIGS)
