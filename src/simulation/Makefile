CXXFLAGS += \
  $(WARNINGS) \
  $(OPTFLAGS) \
  $(DBGFLAGS) \
  $(INCLUDES)

WARNINGS = \
  -std=c++20 \
  -Wpedantic \
  -Wall \
  -Wextra \
  -Wconversion \
  -Wsign-conversion \
  -Wshadow

OPTFLAGS = \
  -funsafe-math-optimizations \
  -O2 \
  -march=native

DBGFLAGS = \
  -g

INCLUDES = \
  -isystem ../../submodules/github.com/snsinfu/cxx-getopt \
  -isystem ../../submodules/github.com/snsinfu/cxx-spline/include \
  -isystem ../../submodules/github.com/snsinfu/h5/include \
  -isystem ../../submodules/github.com/snsinfu/micromd/include \
  -isystem ../../submodules/github.com/snsinfu/tsv/include \
  -isystem ../../submodules/github.com/danielaparker/jsoncons/include

LIBS = \
  -lhdf5

COMMON_SOURCES       = $(wildcard common/*.cpp)
PREPARE_SOURCES      = $(COMMON_SOURCES) $(wildcard stage_prepare/*.cpp)
ANATELOPHASE_SOURCES = $(COMMON_SOURCES) $(wildcard stage_anatelophase/*.cpp)
INTERPHASE_SOURCES   = $(COMMON_SOURCES) $(wildcard stage_interphase/*.cpp)

SOURCES = $(sort $(PREPARE_SOURCES) $(ANATELOPHASE_SOURCES) $(INTERPHASE_SOURCES))
OBJECTS = $(SOURCES:.cpp=.o)

PRODUCTS = \
  prepare \
  anatelophase \
  interphase

ARTIFACTS = \
  $(PRODUCTS) \
  $(OBJECTS) \
  _depends.mk


.PHONY: all clean depends

all: $(PRODUCTS)
	@:

clean:
	rm -f $(ARTIFACTS)

depends:
	for src in $(SOURCES); do \
	    $(CXX) $(CXXFLAGS) -c -MM -MF- -MT $${src%.*}.o $${src}; \
	done > _depends.mk

prepare: $(PREPARE_SOURCES:.cpp=.o)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

anatelophase: $(ANATELOPHASE_SOURCES:.cpp=.o)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

interphase: $(INTERPHASE_SOURCES:.cpp=.o)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

-include _depends.mk
-include _overrides.mk
