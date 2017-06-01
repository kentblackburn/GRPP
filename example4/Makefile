NAME    = example4
CFILES  = allocate.c inverse.c example4.c
CFLAGS  = -O
LDFLAGS = 

LIBS    = -lm

# Rules...

SRCFILES = $(CFILES)
OBJFILES = $(CFILES:.c=.o)

$(NAME): $(OBJFILES) 
	 $(CC) -o $@ $(CFLAGS) $(OBJFILES) $(LIBS) $(LDFLAGS)

example4.c:	example4.g
	grpp -o example4.c example4.g
