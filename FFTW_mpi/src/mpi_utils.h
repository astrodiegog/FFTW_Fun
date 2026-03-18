
/* \fn int greatest_prime_factor(int) */
/* Calculate the greatest prime factor */
extern int greatest_prime_factor(int num);

/* \fn void Tile_Decomposition1D(int, int *) */
/* 1D tile decomposition */
void Tile_Decomposition1D(int nprocs, int *np_x);

void Domain_Decomposition1D(int Nx, int np_x, int *Nx_local);

