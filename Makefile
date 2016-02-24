C=		gcc
CXX=		g++
CFLAGS=		-g -Wall -std=c99 -fopenmp # --stack-check # -fstack-protector-all
CXXFLAGS=	$(CFLAGS)
DFLAGS=		-DHAVE_PTHREAD
DBGFLAGS= 	#-DDEBUG_ENABLED
OBJS=		common.o preprocess_vcf.o #svm/svm.o refine_svm.o subtree.o
PROG=		read_vcf
INCLUDES=	
LIBS=		-lm -lpthread #-lz

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $(DBGFLAGS) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

read_vcf:$(OBJS)
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) -o $@ $(LIBS) #-static
#		$(CXX) $(CFLAGS) $(DFLAGS) $(OBJS) -o $@ $(LIBS) # -static

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a
