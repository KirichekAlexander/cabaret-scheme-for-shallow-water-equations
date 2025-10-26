CXXFLAGS = -std=c++20 -Wall -Wpedantic -fsanitize=undefined,address -Wextra

main: read_matrix.o functions.o

functions.o: functions.h

read_matrix.o: read_matrix.h

clean:
	rm -f main *.o

.PHONY: clean