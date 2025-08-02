CXX       := clang++
CXXFLAGS  := -std=c++17 -I. -Istubs -I../src -pthread
LDFLAGS   := -lgtest -lgtest_main

ifeq ($(COVERAGE),1)
	CXXFLAGS += -fprofile-instr-generate -fcoverage-mapping -DENABLE_COVERAGE
	LDFLAGS  += -fprofile-instr-generate
	PROFILE_OUTPUT := LLVM_PROFILE_FILE=coverage.profraw
endif

SRCS := evaluator_tests.cpp \
		../src/evaluator.cpp \
		../src/nucleotidetree.cpp \
		../src/options.cpp \
		../src/fastareader.cpp \
		../src/read.cpp \
		../src/sequence.cpp \
		../src/knownadapters.cpp \
		stubs/fastqreader.cpp

# Map sources to object files
OBJS := $(patsubst ../src/%.cpp,build/src/%.o,$(SRCS))
OBJS := $(patsubst stubs/%.cpp,build/stubs/%.o,$(OBJS))
OBJS := $(patsubst %.cpp,build/%.o,$(OBJS))

TARGET := evaluator_tests

.PHONY: all clean test coverage-report

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Correctly scoped pattern rules
build/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

build/stubs/%.o: stubs/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

build/src/%.o: ../src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

test: all
	$(PROFILE_OUTPUT) ./$(TARGET)

coverage-report: test
	llvm-profdata merge -sparse coverage.profraw -o coverage.profdata
	llvm-cov show ./$(TARGET) -instr-profile=coverage.profdata \
		-format=html -output-dir=coverage_html \
		-ignore-filename-regex='^/usr'

clean:
	rm -f $(TARGET) coverage.profraw coverage.profdata *.fq
	rm -rf build coverage_html
