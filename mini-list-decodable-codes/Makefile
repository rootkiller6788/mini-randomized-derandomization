CC:=cc
CFLAGS:=-std=c11 -Wall -Wextra -O2 -g -Iinclude -Wno-unused-variable -Wno-unused-function -Wno-unused-parameter
LDFLAGS:=-lm
SRC:=$(wildcard src/*.c)
HEADERS:=$(wildcard include/*.h)

.PHONY: all test check-lines clean examples

all: test examples

test: tests/test.c $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) -o test tests/test.c $(SRC) $(LDFLAGS)
	./test

examples: example1 example2 example3

example1: examples/example1.c $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) -o example1 examples/example1.c $(SRC) $(LDFLAGS)

example2: examples/example2.c $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) -o example2 examples/example2.c $(SRC) $(LDFLAGS)

example3: examples/example3.c $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) -o example3 examples/example3.c $(SRC) $(LDFLAGS)

check-lines:
	@echo "include/: $$(cat include/*.h 2>/dev/null | wc -l) lines"
	@echo "src/:     $$(cat src/*.c 2>/dev/null | wc -l) lines"
	@echo "Total include+src: $$(cat include/*.h src/*.c 2>/dev/null | wc -l) lines"

clean:
	rm -f test example1 example2 example3 *.exe
