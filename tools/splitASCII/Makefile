# Compiler
CXX = mpic++

# Compiler flags
CXXFLAGS = -O3 -std=c++11 


# Source files
SRCS = splitASCII.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable
EXEC = splitASCII

# Include directory
INCLUDES = -I.

.PHONY: all clean

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJS) -o $(EXEC)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)

