include /home/ysw/Downloads/CurveLab-2.1.3/makefile.opt


CC=g++
#CFLAGS+=-g
#CFLAGS+=`pkg-config --cflags opencv`
#LDFLAGS+=`pkg-config --libs opencv`
LDFLAGS+=`pkg-config --libs gtk+-2.0`
CFLAGS+=`pkg-config --cflags gtk+-2.0`


PROG=mainwindow
OBJS=$(PROG).o TV_function.o myfunction.o 

.PHONY: all clean
$(PROG): $(OBJS) libfdct_wrapping.a 
	$(CC) $(CFLAGS) -o $(PROG) $(OBJS)  libfdct_wrapping.a $(LDFLAGS)

%.o: %.cpp
	$(CC) -c $(CFLAGS) $<

all: $(PROG)

clean:
	rm -f $(OBJS) $(PROG)
