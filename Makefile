CXXFLAGS = -std=c++20 -Wall -Wpedantic -fsanitize=undefined,address -Wextra -I./lib/tecio/include
LDFLAGS = -L./lib/tecio -ltecio -fsanitize=undefined,address

main: read_matrix.o functions.o file_manager.o

file_manager.o: file_manager.h

functions.o: functions.h

read_matrix.o: read_matrix.h

clean:
	rm -f main *.o

.PHONY: clean