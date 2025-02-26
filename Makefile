# change this TARGET to compile your own programs
TARGET = main
SHELL  = /bin/sh

DBGFLAGS   = -g -Wall
CFLAGS     = -O3 -Wall -Wno-deprecated -fno-strict-aliasing -fomit-frame-pointer 
EXTRACFLAGS=

CC         = /usr/bin/gcc
CPPFLAGS   = -I. -I/usr/local/include
LD         = /usr/bin/gcc
LDFLAGS    = -O3 -Wall -fno-strict-aliasing -fomit-frame-pointer -Wl,-search_paths_first 
MODLD      = /usr/bin/gcc
MODLDFLAGS = -bundle -undefined dynamic_lookup $(CFLAGS) $(DLCFLAGS)
EXTRAMODLDFLAGS = 
EXTRALIBS  =

RUNPTH     = 
DLCFLAGS   = -fPIC
LIBS       = -lpari -lssl -lcrypto 

RM = rm -f


OBJS = $(TARGET).o
DYN = lib$(TARGET).dylib
ALL = $(TARGET)-sta $(TARGET)-dyn $(DYN)

dft: $(TARGET)-dyn

all: $(ALL)

sta: $(TARGET)-sta

dyn: $(TARGET)-dyn

dynlib: $(DYN)

$(DYN): $(OBJS) $(HEADERS)
	$(MODLD) -o $@ $(MODLDFLAGS) $(EXTRACFLAGS) $(OBJS) $(EXTRAMODLDFLAGS) 

$(TARGET)-sta: $(OBJS)
	$(LD) -o $@ $(LDFLAGS) $(EXTRACFLAGS) $< $(EXTRALIBS) $(STATIC) $(LIBS)

$(TARGET)-dyn: $(OBJS)
	$(LD) -o $@ $(LDFLAGS) $(EXTRACFLAGS) $< $(EXTRALIBS) $(RUNPTH) $(LIBS)

%.o: %.c 
	$(CC) -c $(CFLAGS) $(EXTRACFLAGS) $(CPPFLAGS) $(DLCFLAGS) $<
clean:
	-$(RM) *.o $(ALL)