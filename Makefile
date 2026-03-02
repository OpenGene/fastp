DIR_INC := ./inc
DIR_SRC := ./src
DIR_OBJ := ./obj

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
INCLUDE_DIRS ?=
LIBRARY_DIRS ?=

SRC := $(wildcard ${DIR_SRC}/*.cpp)
OBJ := $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC}))

TARGET := fastp

BIN_TARGET := ${TARGET}

CXX ?= g++
CXXFLAGS := -std=c++11 -pthread -g -O3 -MD -MP -I. -I${DIR_INC} $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir)) ${CXXFLAGS}
LIBS := -lisal -ldeflate -lhwy -lpthread

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
  # Linux: fully static binary
  LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) -static -Wl,--no-as-needed -pthread $(LIBS) $(LD_FLAGS)
else
  # macOS: link 3rd-party libs statically via .a when available, fallback to dynamic
  FIND_STATIC = $(firstword $(foreach d,$(LIBRARY_DIRS),$(wildcard $(d)/lib$(1).a)) $(wildcard /usr/local/lib/lib$(1).a /opt/homebrew/lib/lib$(1).a))
  STATIC_LIBS :=
  DYNAMIC_LIBS :=
  $(foreach lib,isal deflate hwy,\
    $(if $(call FIND_STATIC,$(lib)),\
      $(eval STATIC_LIBS += $(call FIND_STATIC,$(lib))),\
      $(eval DYNAMIC_LIBS += -l$(lib))))
  LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(STATIC_LIBS) $(DYNAMIC_LIBS) -lpthread $(LD_FLAGS)
endif


${BIN_TARGET}:${OBJ}
	$(CXX) $(OBJ) -o $@ $(LD_FLAGS)

${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp
	@mkdir -p $(@D)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

.PHONY:clean
clean:
	@rm -rf $(DIR_OBJ)
	@rm -f $(TARGET)

install:
	install $(TARGET) $(BINDIR)/$(TARGET)
	@echo "Installed."

-include $(OBJ:.o=.d)
