CXX ?= g++

CPPFLAGS ?=
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra -pthread
LDFLAGS ?=
LDLIBS ?= -pthread

SRCS := io_utils.cpp matrix.cpp speedup.cpp fastp.cpp nw_align.cpp main.cpp
OBJS := $(SRCS:.cpp=.o)
TARGET := topn

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
