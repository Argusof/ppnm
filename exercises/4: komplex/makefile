CFLAGS = -Wall
LDLIBS = -lm

.PHONEY: all clean
all              : out.txt           ; cat $<
out.txt          : main              ; ./$< > $@
main             : main.o komplex.o equal.o  # built-in linking rule is used here
main.o komplex.o : komplex.h         # built-in compiling rule is used here
equal.o					 : equal.c funcs.h
clean            :                   ; $(RM) main main.o komplex.o out.txt
