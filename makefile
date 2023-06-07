.PHONY: all

CXX_EXTRA_FLAGS ?= -O3 -march=native -mcpu=native -I./vendored

ifdef DEBUG
CXX_EXTRA_FLAGS += -g3 -DDEBUG
endif
	
VALGRIND_FLAGS ?= --suppressions=valgrind-suppressions.txt  --leak-check=full --show-leak-kinds=all

FFTCONVOLVER_SOURCES = vendored/FFTConvolver/AudioFFT.cpp vendored/FFTConvolver/Utilities.cpp vendored/FFTConvolver/FFTConvolver.cpp

FFTCONVOLVER_OBJECTs = ${FFTCONVOLVER_SOURCES:.cpp=.o}

all: plugins test_eq_match

plugins: lv2/fps-plugins.lv2/dynamics.so lv2/fps-plugins.lv2/eq_match.so

lv2/fps-plugins.lv2/dynamics.so: dynamics.cc makefile
	g++ -std=c++20 ${CXX_EXTRA_FLAGS} dynamics.cc -pedantic -Wall -Werror -shared -o lv2/fps-plugins.lv2/dynamics.so

lv2/fps-plugins.lv2/eq_match.so: eq_match.cc eq_match.h makefile
	g++ -std=c++20 ${CXX_EXTRA_FLAGS} eq_match.cc ${FFTCONVOLVER_SOURCES} -pedantic -Wall -Werror -shared -o lv2/fps-plugins.lv2/eq_match.so -lfftw3f -lm

test_eq_match: test_eq_match.cc eq_match.h makefile
	g++ -std=c++20 ${CXX_EXTRA_FLAGS} test_eq_match.cc ${FFTCONVOLVER_SOURCES} -o test_eq_match -lfftw3f -lm -lsndfile

test: all
	LV2_PATH=${PWD}/lv2 lv2info https://dfdx.eu/fps-plugins.lv2/relative_dynamics
	LV2_PATH=${PWD}/lv2 valgrind ${VALGRIND_FLAGS} lv2bench https://dfdx.eu/fps-plugins.lv2/relative_dynamics
	LV2_PATH=${PWD}/lv2 lv2info https://dfdx.eu/fps-plugins.lv2/eq_match
	LV2_PATH=${PWD}/lv2 valgrind ${VALGRIND_FLAGS} lv2bench https://dfdx.eu/fps-plugins.lv2/eq_match
	valgrind ${VALGRIND_FLAGS} ./test_eq_match 2048 input.wav output.wav linear_phase_response.wav minimum_phase_response.wav matched.wav

