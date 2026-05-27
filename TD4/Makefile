all: parametrization

clean : 
	-rm main.o  parametrization

parametrization: main.o
	g++ -g -o parametrization main.o -O3

main.o : main.cpp
	g++ -c -g main.cpp -O3