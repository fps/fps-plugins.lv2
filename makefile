.PHONY: all

PREFIX ?= /usr/local
INSTALL_DIR ?= ${PREFIX}/lib/lv2

# CXX_EXTRA_FLAGS ?= -O1 -g -fsanitize=address -march=native -mcpu=native -I./vendored -Wall -pedantic
CXX_EXTRA_FLAGS ?= -O3 -march=native -mtune=native -I./vendored -Wall -pedantic `pkg-config lv2 sndfile fftw3f --cflags`
LDFLAGS ?= -fPIC `pkg-config lv2 sndfile fftw3f --libs`

ifdef DEBUG
CXX_EXTRA_FLAGS += -g3 -DDEBUG
endif
	
VALGRIND_FLAGS ?= --suppressions=valgrind-suppressions.txt  --leak-check=full --show-leak-kinds=all

FFTCONVOLVER_SOURCES = vendored/FFTConvolver/AudioFFT.cpp vendored/FFTConvolver/Utilities.cpp vendored/FFTConvolver/FFTConvolver.cpp

all: plugins test_eq_match

plugins: lv2/fps-plugins.lv2/dynamics.so lv2/fps-plugins.lv2/eq_match.so lv2/fps-plugins.lv2/stereo_decorrelation.so

lv2/fps-plugins.lv2/dynamics.so: dynamics.cc makefile
	g++ -std=c++20 ${CXX_EXTRA_FLAGS} dynamics.cc -shared -o lv2/fps-plugins.lv2/dynamics.so ${LDFLAGS}

lv2/fps-plugins.lv2/eq_match.so: eq_match.cc eq_match.h makefile
	g++ -std=c++20 ${CXX_EXTRA_FLAGS} eq_match.cc ${FFTCONVOLVER_SOURCES} -shared -o lv2/fps-plugins.lv2/eq_match.so ${LDFLAGS}

test_eq_match: test_eq_match.cc eq_match.h makefile
	g++ -std=c++20 ${CXX_EXTRA_FLAGS} test_eq_match.cc ${FFTCONVOLVER_SOURCES} -o test_eq_match ${LDFLAGS}

lv2/fps-plugins.lv2/stereo_decorrelation.so: stereo_decorrelation.cc stereo_decorrelation.h makefile
	g++ -std=c++20 ${CXX_EXTRA_FLAGS} stereo_decorrelation.cc ${FFTCONVOLVER_SOURCES} -shared -o lv2/fps-plugins.lv2/stereo_decorrelation.so ${LDFLAGS}
	
test: all
	LV2_PATH=${PWD}/lv2 lv2info https://dfdx.eu/fps-plugins.lv2/relative_dynamics
	LV2_PATH=${PWD}/lv2 valgrind ${VALGRIND_FLAGS} lv2bench https://dfdx.eu/fps-plugins.lv2/relative_dynamics
	LV2_PATH=${PWD}/lv2 lv2info https://dfdx.eu/fps-plugins.lv2/eq_match
	LV2_PATH=${PWD}/lv2 valgrind ${VALGRIND_FLAGS} lv2bench https://dfdx.eu/fps-plugins.lv2/eq_match
	valgrind ${VALGRIND_FLAGS} ./test_eq_match 2048 input.wav output.wav linear_phase_response.wav minimum_phase_response.wav matched.wav

install: all
	mkdir -p ${DESTDIR}/${INSTALL_DIR}
	cp -r lv2/fps-plugins.lv2 ${DESTDIR}/${INSTALL_DIR}
