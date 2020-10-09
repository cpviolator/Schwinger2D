#Your path to

#============================================================
TARGET  = wilson2D
SOURCES = main.cpp
OBJS    = main.o
INC_PATH=-I../include

CXX=g++
CXXFLAGS = -O3 -g -std=c++11 ${INC_PATH}
#============================================================
all: ${TARGET}

${TARGET}: ${OBJS}
	$(CXX) ${CXXFLAGS} -o ${TARGET} ${OBJS} ${ARPACK_LIB}

main.o: main.cpp Makefile 
	${CXX} ${CXXFLAGS} -c main.cpp
#============================================================

clean:
	rm -f ${TARGET} ${OBJS} core*
