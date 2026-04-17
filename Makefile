DIR_INC := ./inc
DIR_SRC := ./src
DIR_OBJ := ./obj

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
INCLUDE_DIRS ?=
LIBRARY_DIRS ?=

# Auto-detect library paths via pkg-config (works with conda, brew, system packages)
PKG_CONFIG ?= pkg-config
HWY_CFLAGS := $(shell $(PKG_CONFIG) --cflags libhwy 2>/dev/null)
HWY_LIBS := $(shell $(PKG_CONFIG) --libs-only-L libhwy 2>/dev/null)
ISAL_CFLAGS := $(shell $(PKG_CONFIG) --cflags libisal 2>/dev/null)
ISAL_LIBS := $(shell $(PKG_CONFIG) --libs-only-L libisal 2>/dev/null)
DEFLATE_CFLAGS := $(shell $(PKG_CONFIG) --cflags libdeflate 2>/dev/null)
DEFLATE_LIBS := $(shell $(PKG_CONFIG) --libs-only-L libdeflate 2>/dev/null)

SRC := $(wildcard ${DIR_SRC}/*.cpp)
OBJ := $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC}))

TARGET := fastp

BIN_TARGET := ${TARGET}

CXX ?= g++
CXXFLAGS := -std=c++11 -pthread -g -O3 -MD -MP -I. -I${DIR_INC} $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir)) $(HWY_CFLAGS) $(ISAL_CFLAGS) $(DEFLATE_CFLAGS) ${CXXFLAGS}
LIBS := -lisal -ldeflate -lhwy -lpthread

PKG_LDFLAGS := $(HWY_LIBS) $(ISAL_LIBS) $(DEFLATE_LIBS)

UNAME_S := $(shell uname -s)
FIND_STATIC = $(firstword $(foreach d,$(LIBRARY_DIRS),$(wildcard $(d)/lib$(1).a)) $(wildcard /usr/local/lib/lib$(1).a /opt/homebrew/lib/lib$(1).a))
STATIC_LIBS :=
DYNAMIC_LIBS :=
$(foreach lib,isal deflate hwy,\
  $(if $(call FIND_STATIC,$(lib)),\
    $(eval STATIC_LIBS += $(call FIND_STATIC,$(lib))),\
    $(eval DYNAMIC_LIBS += -l$(lib))))

ifeq ($(UNAME_S),Linux)
  ifeq ($(DYNAMIC_LIBS),)
    # All .a found: fully static binary (default for Linux)
    LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(PKG_LDFLAGS) -static -Wl,--no-as-needed -pthread $(LIBS) $(LD_FLAGS)
  else
    # Some .a missing (e.g. conda): link .a directly + dynamic fallback
    LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(PKG_LDFLAGS) $(STATIC_LIBS) $(DYNAMIC_LIBS) -lpthread $(LD_FLAGS)
  endif
else
  # macOS: .a preferred, fallback to dynamic
  LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(PKG_LDFLAGS) $(STATIC_LIBS) $(DYNAMIC_LIBS) -lpthread $(LD_FLAGS)
endif


${BIN_TARGET}:${OBJ}
	$(CXX) $(OBJ) -o $@ $(LD_FLAGS)

# to force fully static linking (Linux only)
STATIC_FLAGS := -static -Wl,--no-as-needed -pthread
STATIC_LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(STATIC_FLAGS) $(LIBS) $(STATIC_LD_FLAGS)
static:${OBJ}
	$(CXX) $(OBJ) -o ${BIN_TARGET} $(STATIC_LD_FLAGS)
.PHONY:static

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
