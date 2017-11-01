DIR_INC = ./inc
DIR_SRC = ./src
DIR_OBJ = ./obj
BINDIR=/usr/local/bin

SRC = $(wildcard ${DIR_SRC}/*.cpp)  
OBJ = $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC})) 

TARGET = fastp

BIN_TARGET = ${TARGET}

CC = g++
CFLAGS = -std=c++11 -g -I${DIR_INC}

${BIN_TARGET}:${OBJ}
	$(CC) $(OBJ) -lz -lpthread -o $@
    
${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp make_obj_dir
	$(CC) $(CFLAGS) -O3 -c  $< -o $@
.PHONY:clean
clean:
	rm obj/*.o
	rm $(TARGET)

make_obj_dir:
	@if test ! -d $(DIR_OBJ) ; \
	then \
		mkdir $(DIR_OBJ) ; \
	fi

install:
	install $(TARGET) $(BINDIR)/$(TARGET)
	@echo "Installed."