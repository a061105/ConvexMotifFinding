FLAG = -O3 


all: clean train

train:	OptStep FW PhaseOne ResOut IOandTrans
	g++ $(FLAG) -o train IOandTrans.o OptStep.o Frank_Wolfe.o OptPhaseOne.o ResOut.o  main.cpp
OptStep:
	g++ $(FLAG) -c -o OptStep.o OptStep.cpp
FW:
	g++ $(FLAG) -c -o Frank_Wolfe.o Frank_Wolfe.cpp
PhaseOne:
	g++ $(FLAG) -c -o OptPhaseOne.o OptPhaseOne.cpp
ResOut:
	g++ $(FLAG) -c -o ResOut.o ResOut.cpp
IOandTrans:
	g++ $(FLAG) -c -o IOandTrans.o IOandTrans.cpp


clean:
	rm -f *.o
