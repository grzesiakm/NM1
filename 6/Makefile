CC=gcc
CFLAGS= -O2 -lm -Wall --pedantic
DEPS=
OBJS=main.o
OUT=main

all: $(OUT)

$(OUT): $(OBJS)
	$(CC)  -o $@  $^ $(CFLAGS)
	
.c.o:
	$(CC) -c -o $@ $< $(CFLAGS) 
	
.PHONY: clean	
clean:
	rm -f *o
	rm -f *~
	rm -f $(OUT)
	
run:
	./main

check:
	valgrind ./main
