N ?= 100										# Make default variable

program: main.c
	g++ -o program main.c

test: program							# Have program as prerequisit
	./program $(N)

utest: program
	./program -1

vtest: program
	./program -2
