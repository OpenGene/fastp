DIR_INC := ./inc
DIR_SRC := ./src
DIR_OBJ := ./obj
DIR_HWY := ./third_party/highway

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
INCLUDE_DIRS ?=
LIBRARY_DIRS ?=

SRC := $(wildcard ${DIR_SRC}/*.cpp)
OBJ := $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC}))

# Highway runtime dispatch support
HWY_OBJS := ${DIR_OBJ}/hwy_targets.o ${DIR_OBJ}/hwy_abort.o
OBJ += $(HWY_OBJS)

TARGET := fastp

BIN_TARGET := ${TARGET}

CXX ?= g++
CXXFLAGS := -std=c++11 -pthread -g -O3 -MD -MP -I. -I${DIR_INC} -I${DIR_HWY} $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir)) ${CXXFLAGS}
LIBS := -lisal -ldeflate -lpthread
STATIC_FLAGS := -static -Wl,--no-as-needed -pthread
LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(LIBS) $(LD_FLAGS)
STATIC_LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(STATIC_FLAGS) $(LIBS) $(STATIC_LD_FLAGS)


${BIN_TARGET}:${OBJ}
	$(CXX) $(OBJ) -o $@ $(LD_FLAGS)

static:${OBJ}
	$(CXX) $(OBJ) -o ${BIN_TARGET} $(STATIC_LD_FLAGS)

${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp
	@mkdir -p $(@D)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

# Highway source files for runtime CPU detection and error handling
${DIR_OBJ}/hwy_targets.o:${DIR_HWY}/hwy/targets.cc
	@mkdir -p $(@D)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

${DIR_OBJ}/hwy_abort.o:${DIR_HWY}/hwy/abort.cc
	@mkdir -p $(@D)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

.PHONY:clean
.PHONY:static
clean:
	@rm -rf $(DIR_OBJ)
	@rm -f $(TARGET)

install:
	install $(TARGET) $(BINDIR)/$(TARGET)
	@echo "Installed."

-include $(OBJ:.o=.d)
