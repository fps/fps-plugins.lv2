.PHONY: all

# CXX_EXTRA_FLAGS ?= -O3 -march=native -mcpu=native -g -I .
CXX_EXTRA_FLAGS ?= -O3 -march=native -mcpu=native -g3 -DNDEBUG -I .
VALGRIND_FLAGS ?= --suppressions=valgrind-suppressions.txt  --leak-check=full --show-leak-kinds=all

all: plugins

plugins: lv2/fps-plugins.lv2/dynamics.so lv2/fps-plugins.lv2/eq_match.so

generate:
	lv2-ttl2c -b lv2/fps-plugins.lv2 -o generated 

lv2/fps-plugins.lv2/dynamics.so: dynamics.cc
	g++ -std=c++20 ${CXX_EXTRA_FLAGS} dynamics.cc -pedantic -Wall -Werror -shared -o lv2/fps-plugins.lv2/dynamics.so

lv2/fps-plugins.lv2/eq_match.so: eq_match.cc eq_match.h
	g++ -std=c++20 ${CXX_EXTRA_FLAGS} eq_match.cc FFTConvolver/AudioFFT.cpp FFTConvolver/Utilities.cpp FFTConvolver/FFTConvolver.cpp -pedantic -Wall -Werror -shared -o lv2/fps-plugins.lv2/eq_match.so -lfftw3f -lm

test: plugins
	g++ -std=c++20 ${CXX_EXTRA_FLAGS} test_eq_match.cc -o test_eq_match -lfftw3f -lm
	LV2_PATH=${PWD}/lv2 lv2ls
	LV2_PATH=${PWD}/lv2 lv2info http://dfdx.eu/plugins/relative_dynamics
	LV2_PATH=${PWD}/lv2 valgrind ${VALGRIND_FLAGS} lv2bench http://dfdx.eu/plugins/relative_dynamics
	LV2_PATH=${PWD}/lv2 lv2info http://dfdx.eu/plugins/eq_match
	LV2_PATH=${PWD}/lv2 valgrind ${VALGRIND_FLAGS} lv2bench http://dfdx.eu/plugins/eq_match
	LV2_PATH=${PWD}/lv2 valgrind ${VALGRIND_FLAGS} ./test_eq_match

