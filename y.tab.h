typedef union	{
	 char	letter;
	 char	string[8192];
	} YYSTYPE;
#define	NL	257
#define	SEMI	258
#define	EQ	259
#define	SYMEQ	260
#define	MARKER	261
#define	NUMBER	262
#define	VARIABLE	263
#define	FUNCTION	264
#define	TYPE	265
#define	RANK1	266
#define	RANK2	267
#define	RANK3	268
#define	RANK4	269
#define	COORD	270
#define	INDICE	271
#define	ARRAY	272
#define	VECTOR	273
#define	ANTISYM	274
#define	COMMENT	275
#define	COVDERIV	276
#define	LIEDERIV	277
#define	UMINUS	278


extern YYSTYPE yylval;
