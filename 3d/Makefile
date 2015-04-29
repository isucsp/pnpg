
NAME = $(shell find . -iname '*.c')
NAMES = $(basename $(NAME))

all: $(NAMES)

$(NAMES): *.c
	g++ $@.c -o $@ -pthread -lnetcdf

#g++ -Wall $@.c -o $@ -pthread -lnetcdf

