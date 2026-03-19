#include "mpi_utils.h"

int greatest_prime_factor(int num)
{
    int prime_factor;

    if ((num == 1) || (num == 2))
    {
        return num;
    }

    prime_factor = 2;

    while(1)
    {
        // keep dividing while evenly divisible
		while ( (num % prime_factor) == 0)
		{
			num /= prime_factor;
		}
        
		// cannot divide any further
		if (num == 1)
		{
			break;
		}
        
		// iterate to next prime factor
		prime_factor += 1;
	}
        
	return prime_factor;
}

void Tile_Decomposition1D(int nprocs, int *np_x)
{
	*np_x = nprocs;
	return;
}

extern void Domain_Decomposition1D(int Nx, int np_x, int *Nx_local)
{
	int n;
	n = Nx % np_x;

	/* Only able to evenly split */
	if (!n)
	{
		*Nx_local = Nx / np_x;
	}
	return;
}


