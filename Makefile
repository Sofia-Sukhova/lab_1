.PHONY: all clean

build: make.o func.o
	mkdir build
	mv task build

plot: build
	cp plot.ipynb build
	cp make_task build
	cp time.csv build
	
		
make.o: func.o
	mpicc -o task main.c func.o

func.o: func.c
	gcc -c func.c

clean:
	rm -r build
	rm func.o
