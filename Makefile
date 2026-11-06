.PHONY: all demo test benchmark hardware os os-benchmark
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra
all: softframe/softframe
softframe/softframe: softframe/main.cpp softframe/fast_bvh.h softframe/sobol.h
	$(CXX) $(CXXFLAGS) -march=native -fopenmp $< -o $@
tests/geometry: tests/geometry.cpp softframe/main.cpp softframe/fast_bvh.h softframe/sobol.h
	$(CXX) $(CXXFLAGS) -fopenmp $< -o $@
test: tests/geometry
	./tests/geometry
benchmark: all
	./softframe/softframe
