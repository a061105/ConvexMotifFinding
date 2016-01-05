FLAG = -O3

all: clean train

train:	OptStep FW PhaseOne ResOut 
	g++ $(FLAG) -o train Frank_Wolfe.o OptPhaseOne.o ResOut.o  main.cpp
OptStep:
	g++ $(FLAG) -c -o OptStep.o OptStep.cpp
FW:
	g++ $(FLAG) -c -o Frank_Wolfe.o Frank_Wolfe.cpp
PhaseOne:
	g++ $(FLAG) -c -o OptPhaseOne.o OptPhaseOne.cpp
ResOut:
	g++ $(FLAG) -c -o ResOut.o ResOut.cpp


clean:
	rm -f *.o