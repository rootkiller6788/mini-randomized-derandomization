# Makefile for mini-bpp-rp-zpp-classes
# Builds test suite and examples for BPP/RP/ZPP complexity classes.

CC := cc
CFLAGS := -std=c11 -Wall -Wextra -O2 -g
LDFLAGS := -lm
INCLUDES := -Iinclude

SRC := $(wildcard src/*.c)
HEADERS := $(wildcard include/*.h)
TESTS := $(wildcard tests/*.c)
EXAMPLES := $(wildcard examples/*.c)

TEST_BINS := $(TESTS:.c=)
EXAMPLE_BINS := $(EXAMPLES:.c=)

.PHONY: all test examples clean check-lines

all: test examples

test: $(TEST_BINS)
	@echo "=== Running Tests ==="
	@for t in $(TEST_BINS); do $$t || exit 1; done
	@echo "=== All Tests Passed ==="

examples: $(EXAMPLE_BINS)

%: %.c $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(SRC) $(LDFLAGS)

check-lines:
	@echo "include/ lines: $$(cat include/*.h | wc -l)"
	@echo "src/ lines: $$(cat src/*.c | wc -l)"
	@echo "Total (include + src): $$(cat include/*.h src/*.c | wc -l)"

clean:
	rm -f $(TEST_BINS) $(EXAMPLE_BINS)

# Safety: no filler scan
audit:
	grep -rn "_fn[0-9]" --include="*.c" --include="*.h" . || true
	grep -rn "_aux[0-9]" --include="*.c" --include="*.h" . || true
	grep -rn "TODO\|FIXME\|stub\|placeholder" --include="*.c" --include="*.h" . || true
	@echo "Audit complete"
