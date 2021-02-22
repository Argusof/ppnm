CFLAGS = -O -Wall -std=gnu11
LDLIBS = -lm
OUTPUTFILES = ouput_cmdline.txt error_cmdline.txt
EXECS = cmdline
OBJCTS = cmdline.o
CFILES = cmdline.c
ARGS = 1 2 3 4 5
DISPLAY = cat

default: $(OUTPUTFILES)
	$(DISPLAY) ouput_cmdline.txt
	$(DISPLAY) error_cmdline.txt

$(OUTPUTFILES): $(EXECS)
	./cmdline $(ARGS) 1> ouput_cmdline.txt 2> error_cmdline.txt

$(EXECS): $(OBJECTS)
#	$(CC) $^ -o $@ $(LDLIBS)

# $(OBJECTS): $(CFILES)
# 	$(CC) -c $@



.PHONEY:clean
clean:
	$(RM) input.txt out* error.txt stdio cmdopt cmdline
