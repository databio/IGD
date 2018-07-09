objs = igd_base.o igd_create.o igd_search.o igd.o 
VPATH = src

igd: $(objs)
	cc -o igd $(objs) -lm -lz

igd_base.o: igd_base.c
igd_create.o: igd_create.c
igd_search.o: igd_search.c
igd.o: igd.c

.PHONY: clean
clean:
	rm igd $(objs)
