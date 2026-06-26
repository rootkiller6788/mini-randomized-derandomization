CC:=cc
CFLAGS:=-std=c11 -Wall -Wextra -O2 -g -Iinclude
LDFLAGS:=-lm
SRC:=$(wildcard src/*.c)
all: test
test: tests/test.c $(SRC)
	$(CC) $(CFLAGS) -o test tests/test.c $(SRC) $(LDFLAGS) && ./test
check-lines:
	@echo "Total: $$(cat include/*.h src/*.c 2>/dev/null | wc -l) lines"
clean:
	rm -f test
