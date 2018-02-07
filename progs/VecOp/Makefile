main: main.c utils.c
	g++ -c -o utils.o utils.c
	ar r libutils.a utils.o

	g++ -c -o main.o main.c
	g++ -o prog main.o libutils.a
