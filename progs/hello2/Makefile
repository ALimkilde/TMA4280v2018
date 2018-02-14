main: main.c hello_utils.c goodbye_utils.c
	g++ -c -o hello_utils.o hello_utils.c
	g++ -c -o goodbye_utils.o goodbye_utils.c
	ar r libutils.a hello_utils.o goodbye_utils.o

	g++ -c -o main.o main.c
	g++ -o hello main.o libutils.a
