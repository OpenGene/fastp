DIR_INC := ./inc
DIR_SRC := ./src
DIR_OBJ := ./obj

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
INCLUDE_DIRS ?=
LIBRARY_DIRS ?=

STRIP	?= strip

SRC := $(wildcard ${DIR_SRC}/*.cpp)
OBJ := $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC}))

TARGET := fastp

BIN_TARGET := ${TARGET}

CXX ?= g++
# Optional flags that the user can override by setting CXXFLAGS in the
# env or make argument.  -pthread is a link flag, and serves no purpose
# in the compile command.  It is handled by -lpthread in LIBS.
CXXFLAGS ?= -g -O3 -MD -MP
# Append required flags to standard CXXFLAGS from env
CXXFLAGS += -std=c++11 -I${DIR_INC} $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))
LIBS := -lisal -ldeflate -lpthread
STATIC_FLAGS := -static -Wl,--no-as-needed -pthread
# Append required flags to standard LDFLAGS from env
LDFLAGS += $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(LIBS)
STATIC_LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(STATIC_FLAGS) $(LIBS) $(STATIC_LD_FLAGS)


${BIN_TARGET}:${OBJ}
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)

static:${OBJ}
	$(CXX) $(OBJ) -o ${BIN_TARGET} $(STATIC_LD_FLAGS)

${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp
	@mkdir -p $(@D)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

.PHONY:clean
.PHONY:static
clean:
	@rm -rf $(DIR_OBJ)
	@rm -f $(TARGET)

# Respect DESTDIR for staged installs (used by most package managers).
# DESTDIR is empty by default, so this will install directly to BINDIR
# unless DESTDIR is supplied by the user.
install:
	mkdir -p $(DESTDIR)$(BINDIR)
	install $(TARGET) $(DESTDIR)$(BINDIR)
	@echo "Installed."

# Many package managers use install-strip target if debugging is not enabled
install-strip: install
	$(STRIP) $(DESTDIR)$(BINDIR)/$(TARGET)

-include $(OBJ:.o=.d)
