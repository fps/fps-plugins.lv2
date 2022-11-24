.PHONY: all

all: plugins test

plugins:
	lv2-ttl2c -b lv2/fps-plugins.lv2 -o generated 
	g++ -std=c++20 dynamics.cc -pedantic -Wall -Werror -shared -o lv2/fps-plugins.lv2/dynamics.so
	# gcc eg_amp.c -pedantic -Wall -Werror -shared -o lv2/example.lv2/amp.so
	# gcc eg_exp.c -pedantic -Wall -Werror -shared -o lv2/example.lv2/exp.so

test:
	LV2_PATH=${PWD}/lv2 lv2ls
	LV2_PATH=${PWD}/lv2 lv2info http://fps.io/plugins/relative_dynamics
	LV2_PATH=${PWD}/lv2 valgrind lv2bench http://fps.io/plugins/relative_dynamics
	# LV2_PATH=${PWD}/lv2 valgrind lv2bench http://lv2plug.in/plugins/eg-exp

