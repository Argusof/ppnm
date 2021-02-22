CC = gcc # the C compiler
CFLAGS = -Ofast -fwrapv # options for the C compiler
LDLIBS = -lm # libraries to link

default: out.txt    # default target is to make out.txt
	cat out.txt # show out.txt on the screen

out.txt: epsilon            # out.txt needs epsilon
	./epsilon > out.txt # run epsilon, output goes into out.txt

epsilon: epsilon.o equal.o                  # epsilon needs epsilon.o
	$(CC) epsilon.o equal.o -o epsilon $(LDLIBS) # link epsilon.o into epsilon

epsilon.o: epsilon.c equal.h                 # epsilon.o needs epsilon.c
	$(CC) $(CFLAGS) -c epsilon.c # compile epsilon.c

equal.o: equal.c equal.h
	$(CC) $(CFLAGS) -c equal.c
	
.PHONEY: clean
clean:                              # clean is a phoney target
	$(RM) epsilon.o epsilon equal.o out.txt # clean up the directory

.PHONEY: test
test:                      # test target used for debugging
	echo $(LDLIBS)
	echo $(CC)
	echo $(RM)
