CFLAGS = -O -Wall -std=gnu11
LDLIBS = -lm
OUTPUTFILES = ouput_cmdline.txt error_cmdline.txt ouput_input.txt error_input.txt ouput_fileread.txt error_fileread.txt
INPUTFILES = input.txt
EXECS = cmdline input fileread
OBJCTS = cmdline.o input.o fileread.o
CFILES = cmdline.c input.c fileread.o
ARGS = 1 2 3 4 5
DISPLAY = cat

default: $(OUTPUTFILES)
	$(DISPLAY) ouput_cmdline.txt
	$(DISPLAY) error_cmdline.txt
	$(DISPLAY) ouput_input.txt
	$(DISPLAY) error_input.txt
	$(DISPLAY) ouput_fileread.txt
	$(DISPLAY) error_fileread.txt


$(OUTPUTFILES): $(EXECS) $(INPUTFILES)
	./cmdline $(ARGS) 1> ouput_cmdline.txt 2> error_cmdline.txt
	./input < input.txt 1> ouput_input.txt 2> error_input.txt
	./fileread < input.txt 1> ouput_fileread.txt 2> error_fileread.txt

$(EXECS): $(OBJECTS)
#	$(CC) $^ -o $@ $(LDLIBS)

# $(OBJECTS): $(CFILES)
# 	$(CC) -c $@



.PHONEY:clean
clean:
	$(RM) out* error.txt stdio cmdopt cmdline input fileread
