/* General Relativity PreProcessor */
/*          Version 2.0            */
/*                                 */
/*     Copyright (c) 1992-1998     */
/*     by James Kent Blackburn     */
/*     All Rights Reserved         */
/*                                 */
/* GRPP header file: example1_grpp.h */


#define GRPP_DIMENSION		 2 
#define GRPP_INDICES		 "IiJjKk" 
#define GRPP_COORDINATES	 "XxYy" 


#define COMPONENT		 double 

typedef struct { 
                COMPONENT x,y;
               } TENSOR_c;

typedef struct { 
                COMPONENT X,Y;
               } TENSOR_C;

typedef struct { 
                COMPONENT xx,xy;
                COMPONENT yx,yy;
               } TENSOR_cc;

typedef struct { 
                COMPONENT xX,xY;
                COMPONENT yX,yY;
               } TENSOR_cC;

typedef struct { 
                COMPONENT Xx,Xy;
                COMPONENT Yx,Yy;
               } TENSOR_Cc;

typedef struct { 
                COMPONENT XX,XY;
                COMPONENT YX,YY;
               } TENSOR_CC;

typedef struct { 
                COMPONENT xxx,xxy;
                COMPONENT xyx,xyy;
                COMPONENT yxx,yxy;
                COMPONENT yyx,yyy;
               } TENSOR_ccc;

typedef struct { 
                COMPONENT xxX,xxY;
                COMPONENT xyX,xyY;
                COMPONENT yxX,yxY;
                COMPONENT yyX,yyY;
               } TENSOR_ccC;

typedef struct { 
                COMPONENT xXx,xXy;
                COMPONENT xYx,xYy;
                COMPONENT yXx,yXy;
                COMPONENT yYx,yYy;
               } TENSOR_cCc;

typedef struct { 
                COMPONENT xXX,xXY;
                COMPONENT xYX,xYY;
                COMPONENT yXX,yXY;
                COMPONENT yYX,yYY;
               } TENSOR_cCC;

typedef struct { 
                COMPONENT Xxx,Xxy;
                COMPONENT Xyx,Xyy;
                COMPONENT Yxx,Yxy;
                COMPONENT Yyx,Yyy;
               } TENSOR_Ccc;

typedef struct { 
                COMPONENT XxX,XxY;
                COMPONENT XyX,XyY;
                COMPONENT YxX,YxY;
                COMPONENT YyX,YyY;
               } TENSOR_CcC;

typedef struct { 
                COMPONENT XXx,XXy;
                COMPONENT XYx,XYy;
                COMPONENT YXx,YXy;
                COMPONENT YYx,YYy;
               } TENSOR_CCc;

typedef struct { 
                COMPONENT XXX,XXY;
                COMPONENT XYX,XYY;
                COMPONENT YXX,YXY;
                COMPONENT YYX,YYY;
               } TENSOR_CCC;

typedef struct { 
                COMPONENT xxxx,xxxy;
                COMPONENT xxyx,xxyy;
                COMPONENT xyxx,xyxy;
                COMPONENT xyyx,xyyy;
                COMPONENT yxxx,yxxy;
                COMPONENT yxyx,yxyy;
                COMPONENT yyxx,yyxy;
                COMPONENT yyyx,yyyy;
               } TENSOR_cccc;

typedef struct { 
                COMPONENT xxxX,xxxY;
                COMPONENT xxyX,xxyY;
                COMPONENT xyxX,xyxY;
                COMPONENT xyyX,xyyY;
                COMPONENT yxxX,yxxY;
                COMPONENT yxyX,yxyY;
                COMPONENT yyxX,yyxY;
                COMPONENT yyyX,yyyY;
               } TENSOR_cccC;

typedef struct { 
                COMPONENT xxXx,xxXy;
                COMPONENT xxYx,xxYy;
                COMPONENT xyXx,xyXy;
                COMPONENT xyYx,xyYy;
                COMPONENT yxXx,yxXy;
                COMPONENT yxYx,yxYy;
                COMPONENT yyXx,yyXy;
                COMPONENT yyYx,yyYy;
               } TENSOR_ccCc;

typedef struct { 
                COMPONENT xxXX,xxXY;
                COMPONENT xxYX,xxYY;
                COMPONENT xyXX,xyXY;
                COMPONENT xyYX,xyYY;
                COMPONENT yxXX,yxXY;
                COMPONENT yxYX,yxYY;
                COMPONENT yyXX,yyXY;
                COMPONENT yyYX,yyYY;
               } TENSOR_ccCC;

typedef struct { 
                COMPONENT xXxx,xXxy;
                COMPONENT xXyx,xXyy;
                COMPONENT xYxx,xYxy;
                COMPONENT xYyx,xYyy;
                COMPONENT yXxx,yXxy;
                COMPONENT yXyx,yXyy;
                COMPONENT yYxx,yYxy;
                COMPONENT yYyx,yYyy;
               } TENSOR_cCcc;

typedef struct { 
                COMPONENT xXxX,xXxY;
                COMPONENT xXyX,xXyY;
                COMPONENT xYxX,xYxY;
                COMPONENT xYyX,xYyY;
                COMPONENT yXxX,yXxY;
                COMPONENT yXyX,yXyY;
                COMPONENT yYxX,yYxY;
                COMPONENT yYyX,yYyY;
               } TENSOR_cCcC;

typedef struct { 
                COMPONENT xXXx,xXXy;
                COMPONENT xXYx,xXYy;
                COMPONENT xYXx,xYXy;
                COMPONENT xYYx,xYYy;
                COMPONENT yXXx,yXXy;
                COMPONENT yXYx,yXYy;
                COMPONENT yYXx,yYXy;
                COMPONENT yYYx,yYYy;
               } TENSOR_cCCc;

typedef struct { 
                COMPONENT xXXX,xXXY;
                COMPONENT xXYX,xXYY;
                COMPONENT xYXX,xYXY;
                COMPONENT xYYX,xYYY;
                COMPONENT yXXX,yXXY;
                COMPONENT yXYX,yXYY;
                COMPONENT yYXX,yYXY;
                COMPONENT yYYX,yYYY;
               } TENSOR_cCCC;

typedef struct { 
                COMPONENT Xxxx,Xxxy;
                COMPONENT Xxyx,Xxyy;
                COMPONENT Xyxx,Xyxy;
                COMPONENT Xyyx,Xyyy;
                COMPONENT Yxxx,Yxxy;
                COMPONENT Yxyx,Yxyy;
                COMPONENT Yyxx,Yyxy;
                COMPONENT Yyyx,Yyyy;
               } TENSOR_Cccc;

typedef struct { 
                COMPONENT XxxX,XxxY;
                COMPONENT XxyX,XxyY;
                COMPONENT XyxX,XyxY;
                COMPONENT XyyX,XyyY;
                COMPONENT YxxX,YxxY;
                COMPONENT YxyX,YxyY;
                COMPONENT YyxX,YyxY;
                COMPONENT YyyX,YyyY;
               } TENSOR_CccC;

typedef struct { 
                COMPONENT XxXx,XxXy;
                COMPONENT XxYx,XxYy;
                COMPONENT XyXx,XyXy;
                COMPONENT XyYx,XyYy;
                COMPONENT YxXx,YxXy;
                COMPONENT YxYx,YxYy;
                COMPONENT YyXx,YyXy;
                COMPONENT YyYx,YyYy;
               } TENSOR_CcCc;

typedef struct { 
                COMPONENT XxXX,XxXY;
                COMPONENT XxYX,XxYY;
                COMPONENT XyXX,XyXY;
                COMPONENT XyYX,XyYY;
                COMPONENT YxXX,YxXY;
                COMPONENT YxYX,YxYY;
                COMPONENT YyXX,YyXY;
                COMPONENT YyYX,YyYY;
               } TENSOR_CcCC;

typedef struct { 
                COMPONENT XXxx,XXxy;
                COMPONENT XXyx,XXyy;
                COMPONENT XYxx,XYxy;
                COMPONENT XYyx,XYyy;
                COMPONENT YXxx,YXxy;
                COMPONENT YXyx,YXyy;
                COMPONENT YYxx,YYxy;
                COMPONENT YYyx,YYyy;
               } TENSOR_CCcc;

typedef struct { 
                COMPONENT XXxX,XXxY;
                COMPONENT XXyX,XXyY;
                COMPONENT XYxX,XYxY;
                COMPONENT XYyX,XYyY;
                COMPONENT YXxX,YXxY;
                COMPONENT YXyX,YXyY;
                COMPONENT YYxX,YYxY;
                COMPONENT YYyX,YYyY;
               } TENSOR_CCcC;

typedef struct { 
                COMPONENT XXXx,XXXy;
                COMPONENT XXYx,XXYy;
                COMPONENT XYXx,XYXy;
                COMPONENT XYYx,XYYy;
                COMPONENT YXXx,YXXy;
                COMPONENT YXYx,YXYy;
                COMPONENT YYXx,YYXy;
                COMPONENT YYYx,YYYy;
               } TENSOR_CCCc;

typedef struct { 
                COMPONENT XXXX,XXXY;
                COMPONENT XXYX,XXYY;
                COMPONENT XYXX,XYXY;
                COMPONENT XYYX,XYYY;
                COMPONENT YXXX,YXXY;
                COMPONENT YXYX,YXYY;
                COMPONENT YYXX,YYXY;
                COMPONENT YYYX,YYYY;
               } TENSOR_CCCC;

#undef COMPONENT 
