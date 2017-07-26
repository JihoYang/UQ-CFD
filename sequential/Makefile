#Include files
SOURCES=helper.c init.c boundary_val.c uvp.c main.c visual.c stochastic_galerkin.c

#Compiler
#--------
CC = gcc

CFLAGS = -Wall -Werror -std=gnu99 -g

#Linker flags
#------------
LDFLAGS= -lm

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=sim

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)


$(OBJECTS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@


