.PHONY: all

CXX_EXTRA_FLAGS ?= -O3 -march=native -mcpu=native -g

all: plugins test

plugins: 
	lv2-ttl2c -b lv2/fps-plugins.lv2 -o generated 
	g++ -std=c++20 ${CXX_EXTRA_FLAGS} dynamics.cc -pedantic -Wall -Werror -shared -o lv2/fps-plugins.lv2/dynamics.so

test: plugins
	LV2_PATH=${PWD}/lv2 lv2ls
	LV2_PATH=${PWD}/lv2 lv2info http://fps.io/plugins/relative_dynamics
	LV2_PATH=${PWD}/lv2 valgrind --leak-check=full --show-leak-kinds=all lv2bench http://fps.io/plugins/relative_dynamics

