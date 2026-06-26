CC:=cc
CFLAGS:=-std=c11 -Wall -Wextra -O2 -g -Iinclude
LDFLAGS:=-lm
SRC:=$(wildcard src/*.c)

all: test

test: tests/test.c $(SRC)
	$(CC) $(CFLAGS) -o test tests/test.c $(SRC) $(LDFLAGS) && ./test

check-lines:
	@echo "include/ + src/ = $$(cat include/*.h src/*.c 2>/dev/null | wc -l) lines"
	@echo "Total project = $$(cat include/*.h src/*.c tests/*.c examples/*.c 2>/dev/null | wc -l) lines"

example1: examples/example1.c $(SRC)
	$(CC) $(CFLAGS) -o example1 examples/example1.c $(SRC) $(LDFLAGS) && ./example1

example2: examples/example2.c $(SRC)
	$(CC) $(CFLAGS) -o example2 examples/example2.c $(SRC) $(LDFLAGS) && ./example2

example3: examples/example3.c $(SRC)
	$(CC) $(CFLAGS) -o example3 examples/example3.c $(SRC) $(LDFLAGS) && ./example3

examples: example1 example2 example3

clean:
	rm -f test example1 example2 example3
