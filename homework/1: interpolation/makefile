CFLAGS    +=  $(shell gsl-config --cflags)
LDLIBS    +=  $(shell gsl-config --libs)
CC         = gcc
DISPLAY    = cat
LDLIBS    += -lm
LINOUTPUT  = linOutputDataFile.txt
QUADOUTPUT = quadOutputDataFile.txt
CUBICOUTPUT = cubicOutputDataFile.txt
INPUT      = cos_tabulated_20.txt
EXECS      = main
OBJECTS 	 = main.o linSpline.o binarySearch.o fillArrayFromInputData.o gslIntFunc.o quadSpline.o cubicSpline.o
HDRS    	 = linSpline.h binarySearch.h gslIntFunc.h quadSpline.h cubicSpline.h
FIGS    	 = linPlot.png quadPlot.png cubicPlot.png
TXT        = output.txt

.PHONEY: default
default: $(FIGS)

linPlot.png: linInterpPlot.gpi linOutputDataFile.txt cos_tabulated_20.txt
	pyxplot $<

quadPlot.png: quadInterpPlot.gpi quadOutputDataFile.txt cos_tabulated_20.txt
	pyxplot $<

cubicPlot.png: cubicInterpPlot.gpi cubicOutputDataFile.txt cos_tabulated_20.txt
		pyxplot $<

$(LINOUTPUT) $(QUADOUTPUT) $(CUBICOUTPUT) $(TXT): $(EXECS) $(INPUT)
	./$< $(INPUT) $(LINOUTPUT) $(QUADOUTPUT) $(CUBICOUTPUT) 1> $(TXT)

$(EXECS): $(OBJECTS) $(HDRS)

clean:
	$(RM) $(OUTPUT) $(EXECS) $(OBJECTS) $(FIGS) $(TXT)
