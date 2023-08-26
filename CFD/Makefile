CCXX := g++
IDIR := include
LDIR := lib
ODIR := build
SRCDIR := src
DOCDIR := documents

CXXFLAGS := -I${IDIR}
FFTFLAGS := -I/home/pritam/FFTW/fftw-3.3.8/usr/include
FFTLIB   := -L/home/pritam/FFTW/fftw-3.3.8/usr/lib -lfftw3

TARGET := DendriticGrowth.cpp

_DEPENDENCIES := General.h Fourier_Quadrature.h
DEPENDENCIES  := ${patsubst %,${IDIR}/%,${_DEPENDENCIES}}

_OBJECTFILES := ${patsubst %.h,%.o,${_DEPENDENCIES}} ${patsubst %.cpp,%.o,${TARGET}}
OBJECTFILES  := ${patsubst %,${ODIR}/%,${_OBJECTFILES}}

all: ${OBJECTFILES}
	${CCXX} -std=c++17 -O3 -o run $^

${ODIR}/%.o: %.cpp ${DEPENDENCIES}
	${CCXX} -std=c++17 -O3 -c -o $@ $< ${CXXFLAGS}

${ODIR}/%.o: ${SRCDIR}/%.cpp ${DEPENDENCIES}
	${CCXX} -std=c++17 -O3 -c -o $@ $< ${CXXFLAGS}

${ODIR}/Fourier_Quadrature.o: ${SRCDIR}/Fourier_Quadrature.cpp ${IDIR}/Fourier_Quadrature.h ${IDIR}/General.h
	${CCXX} -std=c++17 -O3 -c -o $@ $< ${FFTFLAGS} ${CXXFLAGS} ${FFTLIB}

run:
	./run 2>/dev/null

.PHONY: clean run document

document:
	@xelatex -synctex=1 -interaction=nonstopmode ${DOCDIR}/Document.tex >> /dev/null
	@rm -rf ${DOCDIR}/*.aux ${DOCDIR}/*.log ${DOCDIR}/*.gz ${DOCDIR}/*.out
	@rm -rf *.aux *.log *.gz *.out
	@mv *.pdf ${DOCDIR}/

clean:
	rm -rf ${ODIR}/*.o *~ run ${IDIR}/*~ Output/*.tec Output/*.vtk Output/*.jpeg
