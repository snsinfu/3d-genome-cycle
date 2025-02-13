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

COMMON_SOURCES       = $(wildcard common/*.cpp) $(wildcard common/forcefield/*.cpp)
PREPARE_SOURCES      = $(COMMON_SOURCES) $(wildcard stage_prepare/*.cpp)
TRANSITION_SOURCES   = $(COMMON_SOURCES) $(wildcard stage_transition/*.cpp)
ANATELOPHASE_SOURCES = $(COMMON_SOURCES) $(wildcard stage_anatelophase/*.cpp)
INTERPHASE_SOURCES   = $(COMMON_SOURCES) $(wildcard stage_interphase/*.cpp)
PROMETAPHASE_SOURCES = $(COMMON_SOURCES) $(wildcard stage_prometaphase/*.cpp)

SOURCES = $(sort \
    $(PREPARE_SOURCES) \
    $(TRANSITION_SOURCES) \
    $(ANATELOPHASE_SOURCES) \
    $(INTERPHASE_SOURCES) \
    $(PROMETAPHASE_SOURCES) \
  )
OBJECTS = $(SOURCES:.cpp=.o)

PRODUCTS = \
  prepare \
  transition \
  anatelophase \
  interphase \
  prometaphase

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

transition: $(TRANSITION_SOURCES:.cpp=.o)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

anatelophase: $(ANATELOPHASE_SOURCES:.cpp=.o)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

interphase: $(INTERPHASE_SOURCES:.cpp=.o)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

prometaphase: $(PROMETAPHASE_SOURCES:.cpp=.o)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

-include _depends.mk
-include _overrides.mk
