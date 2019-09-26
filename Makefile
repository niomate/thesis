# SRC = src
# BIN = bin
# LIB = lib
# INCLUDE = include

# CC=gcc
# CFLAGS=-I$(INCLUDE) -Wall -g -O2 -lm
# CFLAGS_LIB=$(CFLAGS) -fPIC

# BINS = $(filter-out $(BIN)/gui, $(wildcard $(BIN)/*))

# CORNER_TARGETS = $(SRC)/corner_detection.c $(SRC)/mask.c $(SRC)/utils.c

# .PHONY: all lib bin clean

# all: bin lib
# bin: corner_detection inpainting mask
# lib: libmask libinpaint

# inpainting:
# 	$(CC) -o $(BIN)/$@ $(SRC)/$@.c $(CFLAGS)

# corner_detection:
# 	$(CC) -o $(BIN)/$@ $(SRC)/$@.c $(SRC)/utils.c $(SRC)/$@_main.c $(CFLAGS)
	
# mask:
# 	$(CC) -o $(BIN)/$@ $(CORNER_TARGETS) $(SRC)/$@_main.c $(CFLAGS)

# libmask:
# 	$(CC) -shared -o $(LIB)/$@.so $(CORNER_TARGETS) $(CFLAGS_LIB)

# libinpaint:
# 	$(CC) -shared -o $(LIB)/$@.so $(SRC)/inpainting.c $(CFLAGS_LIB)

# clean:
# 	rm $(BINS)
# 	rm $(LIB)/*

CC=gcc
CFLAGS = -Iinclude -Wall -g -O2
LFLAGS = -L/home/danielg/uni/thesis/lib -Wl,-rpath=/home/danielg/uni/thesis/lib
LIBS = -lm -lmask

LIB_TARGETS = $(addprefix src/, corner_detection.c mask_algo.c utils.c)

TARGETS = $(addprefix src/, corner_detection_main.c mask_main.c)

all: lib bin
lib: libmask libinpaint
bin: corner_detection mask inpainting

corner_detection:
	$(CC) $(CFLAGS) $(LFLAGS) -c src/corner_detection_main.c -o bin/corner_detection -lm -lmask

mask:
	$(CC) $(CFLAGS) $(LFLAGS) -c src/mask_main.c -o bin/mask -lm -lmask

inpainting:
	$(CC) $(CFLAGS) -c src/inpainting.c -o bin/inpainting -lm

libmask:
	$(CC) $(CFLAGS) -shared -fPIC $(LIB_TARGETS) -o lib/libmask.so

libinpaint:
	$(CC) $(CFLAGS) -shared -fPIC src/inpainting.c -o lib/libinpaint.so