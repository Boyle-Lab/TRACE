#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

/*
** hmmgetseed() generates an arbitary seed for the random number generator.
*/
int  hmmgetseed(void) 
{
	return ((int) getpid());
}

/* 
** hmmsetseed() sets the seed of the random number generator to a
** specific value.
*/
void hmmsetseed(int seed) 
{
	srand(seed);
}

/*
**  hmmgetrand() returns a (double) pseudo random number in the
**  interval [0,1).
*/

double hmmgetrand(void)
{
	return (double) rand()/RAND_MAX;
}
