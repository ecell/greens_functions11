all: GF2DAbsSym GF3DAbsSym

GF2DAbsSym:
	g++ -std=c++11 -Wall -Wextra -Wpedantic GreensFunction2DAbsSym.cpp -O2 -march=native -mtune=native -I../include/ -I../extlib/ -I../extlib/greens_functions/build ../extlib/greens_functions/build/libgreens_functions.a -lgsl

GF3DAbsSym:
	g++ -std=c++11 -Wall -Wextra -Wpedantic GreensFunction3DAbsSym.cpp -O2 -march=native -mtune=native -I../include/ -I../extlib/ -I../extlib/greens_functions/build ../extlib/greens_functions/build/libgreens_functions.a -lgsl
