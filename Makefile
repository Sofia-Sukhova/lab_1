all: build 
	mkdir build
	mv task build

build:
	mpicc -o task main.c func.c
	
clean:
	rm -r build



