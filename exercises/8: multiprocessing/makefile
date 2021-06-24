CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += -pthread
CFLAGS += -fopenmp
LDLIBS = -lm -lgomp -lpthread

OBJECTS = main.o rndNum.o piApproxMultithreaded.o piApproxOpenMP.o diffClock.o
HDRS = rndNum.h piApproxMultithreaded.h piApproxOpenMP.h diffClock.h
EXECS = main

default: $(EXECS)
	./$<

$(EXECS): $(OBJECTS)

$(OBJECTS): $(HDRS)

clean:
	$(RM) $(EXECS) $(OBJECTS)
