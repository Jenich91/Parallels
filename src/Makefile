.PHONY: all test gcov_report clang cppcheck leaks valgrind clean main winograd gauss ant
CC = g++ -std=c++17
WARNINGS_FLAGS = -Wall -Werror -Wextra
CPPCHECK_FLAGS = --enable=all --language=c++ --std=c++17  --suppress=missingInclude --suppress=missingIncludeSystem --suppress=unusedFunction --suppress=unusedStructMember --suppress=unmatchedSuppression

SOURCES_ALGORITHMS :=winograd/sfleta_*.cc gauss/sfleta_*.cc ant/sfleta_*.cc
HEADERS =:sfleta_*.h winograd/sfleta_*.h gauss/sfleta_*.h ant/sfleta_*.h

OS = $(shell uname)
ifeq ($(OS),Darwin)
	TESTFLAGS = -lgtest
	GCOVFLAGS = --coverage
else
	TESTFLAGS = -lgtest -lpthread -lrt -lm
	GCOVFLAGS = -fprofile-arcs -ftest-coverage
endif

all: clean

clean:
	rm -f *.o *.gcno *.gcda *.info *.out */.clang-format .clang-format
	rm -rf ./report/ *.dSYM
	@clear

gcov_report:
	$(CC) $(WARNINGS_FLAGS) $(GCOVFLAGS) test.cc $(SOURCES_ALGORITHMS) $(TESTFLAGS)  -o gcov_rep.out
	./gcov_rep.out
	lcov -t "gcov_rep" --no-external --output-file=coverag.info --capture --directory .
	genhtml -o ./report coverag.info
	open ./report/index.html

clang:
	@cp ../materials/linters/.clang-format ./
	@clang-format -i *.cc
	@clang-format -i *.h
	clang-format -n *.cc
	clang-format -n *.h

	@cp ../materials/linters/.clang-format ant/
	@clang-format -i ant/*.cc
	@clang-format -i ant/*.h
	clang-format -n ant/*.cc
	clang-format -n ant/*.h
	
	@cp ../materials/linters/.clang-format gauss/
	@clang-format -i gauss/*.cc
	@clang-format -i gauss/*.h
	clang-format -n gauss/*.cc
	clang-format -n gauss/*.h
	
	@cp ../materials/linters/.clang-format winograd/
	@clang-format -i winograd/*.cc
	@clang-format -i winograd/*.h
	clang-format -n winograd/*.cc
	clang-format -n winograd/*.h

cppcheck:
	@cppcheck $(CPPCHECK_FLAGS) $(SOURCES_ALGORITHMS)
	@cppcheck $(CPPCHECK_FLAGS) $(HEADERS)
	@cppcheck $(CPPCHECK_FLAGS) sfleta_interface.cc
	@cppcheck $(CPPCHECK_FLAGS) main.cc

leaks:
	$(CC) $(WARNINGS_FLAGS) $(TESTFLAGS) $(SOURCES_ALGORITHMS) test.cc -o test.out
	CK_FORK=no leaks --atExit -- ./test.out

valgrind:
	valgrind --leak-check=full --show-leak-kinds=all --trace-children=yes ./test.out

test: clean
	$(CC) $(WARNINGS_FLAGS) test.cc $(SOURCES_ALGORITHMS) $(TESTFLAGS) -o test.out
	./test.out

main:
	@clear
	$(CC) $(WARNINGS_FLAGS) $(SOURCES_ALGORITHMS) sfleta_interface.cc main.cc -o main.out

winograd: main
	./main.out "winograd"

gauss: main
	./main.out "gauss"

ant: main
	./main.out "ant"

