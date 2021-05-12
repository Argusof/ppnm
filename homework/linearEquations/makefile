CFLAGS    +=  $(shell gsl-config --cflags)
LDLIBS    +=  $(shell gsl-config --libs)
CC         = gcc
DISPLAY    = cat
LDLIBS    += -lm
OUTPUT     = outputFile.txt gslOutput.txt
EXECS      = main
OBJECTS 	 = main.o GSdecomp.o backsub.o utilities.o diffClock.o test_runtime.o
HDRS    	 = GSdecomp.h backsub.h utilities.h diffClock.h
FIGS    	 = runTimePlot.png runTimePlotGSL.png
TXT        = output.txt

default: $(TXT) $(FIGS)
	cat $<

$(TXT): $(EXECS)
	./$< $(OUTPUT) 1> $@

runTimePlot.png: runTimePlot.gpi outputFile.txt
	pyxplot $<

runTimePlotGSL.png: runTimePlotGSL.gpi gslOutput.txt
	pyxplot $<

$(EXECS): $(OBJECTS) $(HDRS)

clean:
	$(RM) $(EXECS) $(OBJECTS) $(TXT) $(OUTPUT) $(FIGS)
