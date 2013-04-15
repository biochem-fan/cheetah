/*
 * Algorithm from N. Wirth's book, implementation by N. Devillard.
 * This code in public domain.
 */

/*
 *	Find kth smallest element of a data array
 *	See: http://ndevilla.free.fr/median/median.pdf
 *	http://ndevilla.free.fr/median/median/src/wirth.c
 *	Reference: N. Wirth: "Algorithms + data structures = programs" Englewood Cliffs: Prentice-Hall, 1976, p. 366
 */

#define median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))

int16_t kth_smallest(int16_t *a, long n, long k);
float kth_smallest(float *a, long n, long k);
