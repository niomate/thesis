CC = gcc
CC_FLAGS = -Wall -Wextra -O3 -g -lm
CC_OBJ_FLAGS = $(CC_FLAGS) -c

SRC_PATH = src
OBJ_PATH = obj

CLEAN_LIST = corner inpaint $(OBJ_PATH) 

.PHONY: all clean

default: all
all: corner inpaint

corner: $(OBJ_PATH)/corner_detection.o
	$(CC) $(CC_FLAGS) -o $@ $?

inpaint: $(OBJ_PATH)/inpainting.o
	$(CC) $(CC_FLAGS) -o $@ $?

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.c
	if [ ! -d $(OBJ_PATH) ]; then\
 	    echo "Create directory ${OBJ_PATH}";\
 	    mkdir -p $(OBJ_PATH);\
	fi
	$(CC) $(CC_OBJ_FLAGS) -o $@ $<

clean:
	@echo CLEAN $(CLEAN_LIST)
	@rm -rf $(CLEAN_LIST)

