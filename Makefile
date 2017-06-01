#	File: 		Makefile
#	Description:	Make file for GRPP
#	Author:		James Kent Blackburn
#	Date:		May 1999

#	Copyright (c) 1992-1999 
#	By James Kent Blackburn
#	All Rights Reserved

CC = 		gcc
CFLAGS = 	-g
LEX =		flex
LFLAGS =	
YACC =		bison
YFLAGS = 	-d -y 
#LDFLAGS = 	-lfl
LDFLAGS =	 
OBJFILES = 	grpp.o grpp_lex.o grpp_yacc.o

install:	grpp clean

test:	grpp
	grpp -o sample.c sample.g

grpp:	$(OBJFILES) 
	$(CC) $(CFLAGS) $(OBJFILES) $(LDFLAGS) -o $@

grpp.o:	grpp.c
	$(CC) $(CFLAGS) -c $<

grpp_lex.o:	grpp_lex.l y.tab.h
	$(LEX) $(LFLAGS) grpp_lex.l
	mv lex.yy.c grpp_lex.c
	$(CC) $(CFLAGS) -c grpp_lex.c

grpp_yacc.o y.tab.h:	grpp_yacc.y
	$(YACC) $(YFLAGS) grpp_yacc.y
	mv y.tab.c grpp_yacc.c
	$(CC) $(CFLAGS) -c grpp_yacc.c

clean:
	rm grpp_yacc.o
	rm grpp_lex.o
	rm grpp.o

realclean:
	rm y.tab.h
	rm grpp_lex.c
	rm grpp_yacc.c
