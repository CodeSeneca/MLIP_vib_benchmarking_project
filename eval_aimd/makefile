.PHONY: all run checkstyle clean compile

TARGET = eval_aimd

CXX = g++
CXXFLAGS = -pedantic -Wall -Wextra -Wshadow -g
INCLUDE_PATHS =
LIBRARY_PATHS =
LDFLAGS =

CPPFILES = $(wildcard *.cpp)
HEADERFILES = $(wildcard *.h)
CPPO = $(patsubst %.cpp, %.o, $(CPPFILES))

all: compile checkstyle run

run:
	./$(TARGET)

checkstyle:
	cpplint --repository=. *.h *.cpp

clean:
	rm -f *.o
	rm -f *.out

compile: $(TARGET)

$(TARGET): $(CPPO)
	$(CXX) $(CXXFLAGS) $(LIBRARY_PATHS) $^ $(LDFLAGS) -o $@

%.o: %.cpp $(HEADERFILES)
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATHS) -c $< -o $@
