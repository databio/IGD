# CFLAGS = -g -lz -lm -O2  # for production/benchmarking
CFLAGS = -g -lz -lm  # for debug
BIN = bin
OBJ = obj
VPATH = src src/libigd
LIB = igd_base.o igd_create.o igd_search.o igd.o 
OBJS = $(addprefix $(OBJ)/, $(LIB))

$(OBJ)/%.o: %.c
	cc -c $(CFLAGS) $< -o $@ 

igd_dev1: $(OBJS)
	cc -o $(BIN)/igd $(OBJS) $(CFLAGS)
all: $(OBJS)

$(OBJS): | $(OBJ) $(BIN)

$(OBJ):
	mkdir -p $(OBJ)

$(BIN):
	mkdir -p $(BIN)

.PHONY: clean
clean:
	rm -rf $(BIN)/*
	rm -rf $(OBJ)/*
	rm -rf src/libigd/*.o
