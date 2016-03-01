run:
	mpicc fdmGnu.c -D GNUPLOT_ONE
	mpirun -n 2 ./a.out 16 30 1000
	gnuplot gnu -p
