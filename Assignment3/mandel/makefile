all: mandel_seq mandel_par

mandel_seq: mandel_seq.c pngwriter.c
	gcc -o $@ -I. -O3 $^ -lpng -lm 

mandel_par: mandel_par.c pngwriter.c
	gcc -o $@ -I. -O3 $^ -lpng -lm -fopenmp


clean:
	rm -rf mandel_seq mandel_par

