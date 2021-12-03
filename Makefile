pimachine: src/pimachine.cc
	g++ -Wall -Wextra src/pimachine.cc -lgmp -o pimachine -O2

debug: src/pimachine.cc
	g++ -Wall -Wextra src/pimachine.cc -g -lgmp -o pimachine -O0
