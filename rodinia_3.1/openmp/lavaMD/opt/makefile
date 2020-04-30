# Example
# target: dependencies
	# command 1
	# command 2
          # .
          # .
          # .
	# command n

CC = gcc
C_C_1 = -g -O3 -Ofast -march=native -mtune=native -fprofile-generate
C_C_2 = -g -O3 -Ofast -march=native -mtune=native -fprofile-use -fprofile-correction
OMP_LIB = -lgomp
OMP_FLAG = -fopenmp -g


all: first run second

first: main1 kernel_cpu1 num1 timer1 lavaMD_1

second: main2 kernel_cpu2 num2 timer2 lavaMD_2


# link objects (binaries) together
lavaMD_1:		main.o \
			./kernel/kernel_cpu.o \
			./util/num/num.o \
			./util/timer/timer.o
	$(CC) $(C_C_1)	main.o \
			./kernel/kernel_cpu.o \
			./util/num/num.o \
			./util/timer/timer.o \
			-lm \
			$(OMP_LIB) \
			-o lavaMD
	./lavaMD -cores 1 -boxes1d 10

lavaMD_2:		main.o \
			./kernel/kernel_cpu.o \
			./util/num/num.o \
			./util/timer/timer.o
	$(CC) $(C_C_2)	main.o \
			./kernel/kernel_cpu.o \
			./util/num/num.o \
			./util/timer/timer.o \
			-lm \
			$(OMP_LIB) \
			-o lavaMD


# compile function files into objects (binaries)
main1:		main.h \
			main.c \
			./kernel/kernel_cpu.h \
			./kernel/kernel_cpu.c \
			./util/num/num.h \
			./util/num/num.c \
			./util/timer/timer.h \
			./util/timer/timer.c
	$(CC) $(C_C_1)	main.c \
			-c \
			-o main.o \
			-O3

main2:		main.h \
			main.c \
			./kernel/kernel_cpu.h \
			./kernel/kernel_cpu.c \
			./util/num/num.h \
			./util/num/num.c \
			./util/timer/timer.h \
			./util/timer/timer.c
	$(CC) $(C_C_2)	main.c \
			-c \
			-o main.o \
			-O3

kernel_cpu1:	./kernel/kernel_cpu.h \
						./kernel/kernel_cpu.c
	$(CC) $(C_C_1)				./kernel/kernel_cpu.c \
						-c \
						-o ./kernel/kernel_cpu.o \
						-O3 \
						$(OMP_FLAG)

kernel_cpu2:	./kernel/kernel_cpu.h \
						./kernel/kernel_cpu.c
	$(CC) $(C_C_2)				./kernel/kernel_cpu.c \
						-c \
						-o ./kernel/kernel_cpu.o \
						-O3 \
						$(OMP_FLAG)

num1:	./util/num/num.h \
					./util/num/num.c
	$(CC) $(C_C_1)			./util/num/num.c \
					-c \
					-o ./util/num/num.o \
					-O3

num2:	./util/num/num.h \
					./util/num/num.c
	$(CC) $(C_C_2)			./util/num/num.c \
					-c \
					-o ./util/num/num.o \
					-O3

timer1:	./util/timer/timer.h \
						./util/timer/timer.c
	$(CC) $(C_C_1)				./util/timer/timer.c \
						-c \
						-o ./util/timer/timer.o \
						-O3

timer2:	./util/timer/timer.h \
						./util/timer/timer.c
	$(CC) $(C_C_2)				./util/timer/timer.c \
						-c \
						-o ./util/timer/timer.o \
						-O3
# delete all object and executable files
clean:
	rm -rf	*.o \
		./kernel/*.o \
		./util/num/*.o \
		./util/timer/*.o \
                ./util/device/*.o \
		main.gcda \
		lavaMD.gcda \
		lavaMD
