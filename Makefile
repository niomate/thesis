CC = gcc
CC_FLAGS = -Wall -Wextra -O3 -g -lm
CC_OBJ_FLAGS = $(CC_FLAGS) -c

SRC_PATH = src
OBJ_PATH = obj
TARGET_PATH = bin

CLEAN_LIST = $(TARGET_PATH) \
             $(OBJ_PATH) 

default: all
.PHONY: all clean

$(TARGET_PATH)/corner: $(OBJ_PATH)/corner_detection.o $(OBJ_PATH)/utils.o 
	if [ ! -d $(TARGET_PATH) ]; then\
		echo "Create directory ${TARGET_PATH}";\
		mkdir -p $(TARGET_PATH);\
	fi
	$(CC) $(CC_FLAGS) -o $@ $?

$(TARGET_PATH)/inpaint: $(OBJ_PATH)/inpainting.o
	if [ ! -d $(TARGET_PATH) ]; then\
		echo "Create directory ${TARGET_PATH}";\
		mkdir -p $(TARGET_PATH);\
	fi
	$(CC) $(CC_FLAGS) -o $@ $?

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.c
	if [ ! -d $(OBJ_PATH) ]; then\
 	    echo "Create directory ${OBJ_PATH}";\
 	    mkdir -p $(OBJ_PATH);\
	fi
	$(CC) $(CC_OBJ_FLAGS) -o $@ $<

all: $(TARGET_PATH)/corner $(TARGET_PATH)/inpaint

clean:
	@echo CLEAN $(CLEAN_LIST)
	@rm -rf $(CLEAN_LIST)

