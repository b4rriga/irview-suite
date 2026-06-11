all:
	gcc -O3 -march=native -ffast-math -shared -fPIC -o ircore.so ircore.c -lm
