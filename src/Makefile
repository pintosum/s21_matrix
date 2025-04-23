.PHONY: all clean gcov_report memcheck test_build check-format format

CC=gcc
CFLAGS=-std=c11 -Wall -Wextra -Werror
CFLAGS_GCOV=$(CFLAGS) -fprofile-arcs -ftest-coverage

SRC=s21_matrix.c
TEST_SRC=test.c
TEST_EXE=test.out

LIB=s21_matrix.a
LIB_OBJ=s21_matrix.o

LDFLAGS := -lcheck -lm -lsubunit
LDFLAGS_GCOV := $(LDFLAGS) -fprofile-arcs --coverage
LIBLINK := -L. -l:$(LIB)
MEMCHECK := valgrind --tool=memcheck --leak-check=yes ./${TEST_EXE}



UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
	LDFLAGS := -lcheck
	LDFLAGS_GCOV := $(LDFLAGS) -fprofile-arcs --coverage
	LIBLINK := $(LIB)
	MEMCHECK := leaks -atExit -- ./${TEST_EXE} | grep LEAK:
endif

all: clean $(LIB)

$(LIB): clean $(LIB_OBJ)
	ar rcs $@ $(LIB_OBJ)
	rm -rf $(LIB_OBJ)

$(LIB_OBJ): clean
	$(CC) $(CFLAGS) $(SRC) -c -o $(LIB_OBJ)

test_build: $(LIB)
	$(CC) $(CFLAGS) $(TEST_SRC) -o $(TEST_EXE) $(LDFLAGS) $(LIBLINK)

test: test_build
	./${TEST_EXE}

gcov_report: clean
	$(CC) $(CFLAGS_GCOV) $(SRC) -c -o $(LIB_OBJ)
	ar rcs $(LIB) $(LIB_OBJ)
	$(CC) $(CFLAGS) $(TEST_SRC) -o $(TEST_EXE) $(LDFLAGS_GCOV) $(LIBLINK)
	./${TEST_EXE}
	lcov -t "s21_matrix" -o coverage.info -c -d .
	genhtml -o report coverage.info

	@echo "\033[32m\t\t                      GCOV REPORT IS DONE                        \033[0m"
	@echo "\033[32m\t\t          report generated in ./report/index.html                \033[0m"
	rm -f *.gcda *.gcno *.o coverage.info


clean:
	rm -rf coverage.info report/ *.o *.a *.out *.gcda *.gcno *.gcov

memcheck: test_build
	${MEMCHECK}

format:
	clang-format -i ${SRC} ${TEST_SRC} *.h

check-format:
	clang-format -n ${SRC} ${TEST_SRC} *.h
