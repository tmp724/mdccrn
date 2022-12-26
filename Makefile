#COMPILER := clang++
COMPILER := g++

COMPILEFLAGS := -std=c++11 -Wall --pedantic

LINKERFLAGS := -lginac -lcln

OBJECTS := src/command_line/crn.o src/command_line/monotonedependenciescalculator.o src/command_line/main.o

.PHONY: all

all: command_line

command_line: src/command_line/math.hpp crn.o monotonedependenciescalculator.o main.o
	$(COMPILER) $(COMPILEFLAGS) src/command_line/math.hpp src/command_line/crn.o src/command_line/monotonedependenciescalculator.o src/command_line/main.o $(LINKERFLAGS) -o bin/linux/mdccrn_command_line

main.o: src/command_line/math.hpp src/command_line/main.cpp
	$(COMPILER) -c $(COMPILEFLAGS) src/command_line/main.cpp -o src/command_line/main.o

monotonedependenciescalculator.o: src/command_line/monotonedependenciescalculator.cpp src/command_line/monotonedependenciescalculator.hpp src/command_line/math.hpp
	$(COMPILER) -c $(COMPILEFLAGS) src/command_line/monotonedependenciescalculator.cpp -o src/command_line/monotonedependenciescalculator.o

crn.o: src/command_line/crn.cpp src/command_line/crn.hpp src/command_line/math.hpp
	$(COMPILER) -c $(COMPILEFLAGS) src/command_line/crn.cpp -o src/command_line/crn.o

clean:
	rm bin/linux/* $(OBJECTS)
