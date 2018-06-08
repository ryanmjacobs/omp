HEADERS=ompDeblur.h  seqDeblur.h

ifdef GPROF
CFLAGS += -O2 -pg
GPROFFLAG = -pg
else
CFLAGS += -O3
endif

CFLAGS += -fopenmp -I. -std=c89
#CFLAGS += -I.
ifndef $(SRCFILE)
SRCFILE=ompDeblur.c	
endif

EXEC_FILE=deblurTest
OBJ_FILE=$(patsubst %.c,%.o,$(SRCFILE))

all: u.file f.file $(EXEC_FILE)

u.file:
	gzip -kfd u.file.gz

f.file:
	gzip -kfd f.file.gz

$(EXEC_FILE): $(OBJ_FILE) main.o seqDeblur.o
	gcc -o $(EXEC_FILE) -fopenmp $(OBJ_FILE) main.o seqDeblur.o -lm $(GPROFFLAG) 

$(OBJ_FILE): $(SRCFILE) $(HEADERS)
	gcc $(CFLAGS) -c $(SRCFILE)

main.o: main.c $(HEADERS)
	gcc $(CFLAGS) -c main.c

seqDeblur.o: seqDeblur.c $(HEADERS)
	gcc $(CFLAGS) -c seqDeblur.c

clean:
	rm -f {u,f}.file
	rm -f $(EXEC_FILE) $(OBJ_FILE) main.o seqDeblur.o gmon.out
