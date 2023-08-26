CCXX := g++
IDIR := include
LDIR := lib
ODIR := build
SRCDIR := src
DOCDIR := documents

CXXFLAGS := -I${IDIR}

TARGET := r-Adaption.cpp

_DEPENDENCIES := General.h
DEPENDENCIES  := ${patsubst %,${IDIR}/%,${_DEPENDENCIES}}

_OBJECTFILES := ${patsubst %.h,%.o,${_DEPENDENCIES}} ${patsubst %.cpp,%.o,${TARGET}}
OBJECTFILES  := ${patsubst %,${ODIR}/%,${_OBJECTFILES}}

all: ${OBJECTFILES}
	${CCXX} -std=c++17 -O3 -o run $^

${ODIR}/%.o: %.cpp ${DEPENDENCIES}
	${CCXX} -std=c++17 -O3 -c -o $@ $< ${CXXFLAGS}

${ODIR}/%.o: ${SRCDIR}/%.cpp ${DEPENDENCIES}
	${CCXX} -std=c++17 -O3 -c -o $@ $< ${CXXFLAGS}

run:
	./run 2>/dev/null

.PHONY: clean run document

document:
	@xelatex -synctex=1 -interaction=nonstopmode ${DOCDIR}/Document.tex >> /dev/null
	@rm -rf ${DOCDIR}/*.aux ${DOCDIR}/*.log ${DOCDIR}/*.gz ${DOCDIR}/*.out
	@rm -rf *.aux *.log *.gz *.out
	@mv *.pdf ${DOCDIR}/

clean:
	rm -rf ${ODIR}/*.o *~ run ${IDIR}/*~
	rm -rf Output/Flame/*.tec Output/Refinement/*.tec Output/Redistribution/*.tec Output/p-Adaption/*.tec
	rm -rf Output/Flame/*.jpeg Output/Refinement/*.jpeg Output/Redistribution/*.jpeg Output/p-Adaption/*.jpeg	
