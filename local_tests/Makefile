CXX       := clang++
CXXFLAGS  := -std=c++17 -I. -Istubs -I../src -pthread
LDFLAGS   := -lgtest -lgtest_main

SRCS      := evaluator_tests.cpp \
             ../src/evaluator.cpp \
             ../src/nucleotidetree.cpp \
             ../src/options.cpp \
             ../src/fastareader.cpp \
             stubs/fastqreader.cpp

# Convert source path to build path — e.g., ../src/foo.cpp → build/../src/foo.o
OBJS      := $(SRCS:.cpp=.o)
OBJS      := $(addprefix build/,$(OBJS))

TARGET    := evaluator_tests

.PHONY: all clean test

all: build $(TARGET)

# Compile each source to a namespaced object file in build/
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

build/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

build/../src/%.o: ../src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

build/stubs/%.o: stubs/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

build:
	mkdir -p build

test: all
	./$(TARGET)

clean:
	rm -f $(TARGET) build/**/*.o build/*.o *.fq
