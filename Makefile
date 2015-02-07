CFLAGS = -g -c -Wall -pedantic -Dlinux --std=c99 -D_POSIX_C_SOURCE=200112L -D_BSD_SOURCE
LFLAGS = -lcpgplot -lpgplot -lm -lgsl

# Mac OS X (with gcc, PGPLOT installed via fink)
ifeq ($(shell uname),Darwin)
    LINKER = gfortran
	CFLAGS += -D_DARWIN_C_SOURCE -I/usr/local/include/
    LFLAGS += -L/usr/X11R6/lib -lX11 -Wl,-framework -Wl,Foundation -lpng -L/usr/local/lib
else
	CC = gcc
	LINKER = gcc
endif

SRC = microlensing.c searchgrid.c typedefs.c
OBJ = $(SRC:.c=.o)

raytrace: $(OBJ)
	$(LINKER) -o $@ $(OBJ) $(LFLAGS)

clean:
	-rm $(OBJ) raytrace

.SUFFIXES: .c
.c.o:
	$(CC) $(CFLAGS) -c $< -o $@