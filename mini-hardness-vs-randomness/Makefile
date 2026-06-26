CC:=cc
CFLAGS:=-std=c11 -Wall -Wextra -O2 -g -Iinclude
LDFLAGS:=-lm
SRC:=$(wildcard src/*.c)
EXAMPLES:=$(wildcard examples/*.c)
BENCHES:=$(wildcard benches/*.c)

all: test examples

test: tests/test.c $(SRC)
	$(CC) $(CFLAGS) -o test tests/test.c $(SRC) $(LDFLAGS) && ./test

examples: $(EXAMPLES)
	@for f in $(EXAMPLES); do \
		name=$$(basename $$f .c); \
		echo "Building $$name..."; \
		$(CC) $(CFLAGS) -o $$name $$f $(filter-out %/test.c,$(SRC)) $(LDFLAGS) && ./$$name; \
	done

benches: $(BENCHES)
	@for f in $(BENCHES); do \
		name=$$(basename $$f .c); \
		echo "Bench: $$name..."; \
		$(CC) $(CFLAGS) -o $$name $$f $(filter-out %/test.c,$(SRC)) $(LDFLAGS) && ./$$name; \
	done

check-lines:
	@echo "Total include/: $$(cat include/*.h 2>/dev/null | wc -l) lines"
	@echo "Total src/:     $$(cat src/*.c src/*.lean 2>/dev/null | wc -l) lines"
	@echo "Combined:       $$(cat include/*.h src/*.c src/*.lean 2>/dev/null | wc -l) lines"

clean:
	rm -f test
	@for f in $(EXAMPLES); do name=$$(basename $$f .c); rm -f $$name $$name.exe; done
	@for f in $(BENCHES); do name=$$(basename $$f .c); rm -f $$name $$name.exe; done
