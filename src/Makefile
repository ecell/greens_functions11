all: GF2DAbsSym GF3DAbsSym GF3DRadInf GF3DRadAbs GF2DRadAbs

GF2DAbsSym: GreensFunction2DAbsSym.cpp ../include/gf11/*
	g++ -std=c++11 -Wall -Wextra -Wpedantic GreensFunction2DAbsSym.cpp -O2 -march=native -mtune=native -I../include/ -I../extlib/ -I../extlib/greens_functions/build ../extlib/greens_functions/build/libgreens_functions.a -lgsl -o GF2DAbsSym

GF3DAbsSym: GreensFunction3DAbsSym.cpp ../include/gf11/*
	g++ -std=c++11 -Wall -Wextra -Wpedantic GreensFunction3DAbsSym.cpp -O2 -march=native -mtune=native -I../include/ -I../extlib/ -I../extlib/greens_functions/build ../extlib/greens_functions/build/libgreens_functions.a -lgsl -o GF3DAbsSym

GF3DRadInf: GreensFunction3DRadInf.cpp ../include/gf11/*
	g++ -std=c++11 -Wall -Wextra -Wpedantic GreensFunction3DRadInf.cpp -O2 -march=native -mtune=native -I../include/ -I../extlib/ -I../extlib/greens_functions/build ../extlib/greens_functions/build/libgreens_functions.a -lgsl -o GF3DRadInf

GF3DRadAbs: GreensFunction3DRadAbs.cpp ../include/gf11/*
	g++ -std=c++11 -Wall -Wextra -Wpedantic GreensFunction3DRadAbs.cpp -O2 -march=native -mtune=native -I../include/ -I../extlib/ -I../extlib/greens_functions/build ../extlib/greens_functions/build/libgreens_functions.a -lgsl -o GF3DRadAbs

GF2DRadAbs: GreensFunction2DRadAbs.cpp ../include/gf11/*
	g++ -std=c++11 -Wall -Wextra -Wpedantic GreensFunction2DRadAbs.cpp -Og -g -march=native -mtune=native -I../include/ -I../extlib/ -I../extlib/greens_functions/build ../extlib/greens_functions/build/libgreens_functions.a -lgsl -o GF2DRadAbs

.PHONY: clean
clean:
	rm GF3DRadInf GF3DAbsSym GF2DAbsSym GF3DRadAbs GF2DRadAbs
