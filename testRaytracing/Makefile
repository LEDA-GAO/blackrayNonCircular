CXX = clang++
#CXXFLAGS = -std=c++17 -Wall -O3 -Iinclude
CXXFLAGS = -std=c++17 -Wall -g -O0 -Iinclude


EXE = testRayTracing

SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:.cpp=.o)

all: $(EXE)

$(EXE): $(OBJ)
	$(CXX) $^ -o $@

%.o: %.cpp include/def.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXE)