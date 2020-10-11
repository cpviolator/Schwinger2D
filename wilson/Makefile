TARGET  = wilson2D
SOURCES = wilson2D.cpp
OBJS    = wilson2D.o 
INC_PATH=-I../include

CXX=g++
CXXFLAGS = -O3 -g -std=c++11 ${INC_PATH} ${ARPACK_FLAGS}
#============================================================

all: ${TARGET}

${TARGET}: ${OBJS}
	$(CXX) ${CXXFLAGS} -o ${TARGET} ${OBJS} ${ARPACK_LIB}

${TARGET}.o: ${TARGET}.cpp Makefile 
	${CXX} ${CXXFLAGS} -c ${TARGET}.cpp
#============================================================

clean:
	rm -f ${TARGET} ${OBJS} core*
