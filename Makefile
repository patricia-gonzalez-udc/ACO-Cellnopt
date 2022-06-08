
all: run_tests

build:
	mkdir build

build/Makefile: build 
	cd build && cmake ../

build/%: build/Makefile
	cd build && make $*

aco_toymodel: build/aco benchmark/toymodel/toymodel.h5
	$^

clean:
	rm -fr build/
