DEBUG=1 PATH=../lv2-ttl2c/:$PATH LD_LIBRARY_PATH=$(dirname $(which lv2ls))/../lib make -j4 test
