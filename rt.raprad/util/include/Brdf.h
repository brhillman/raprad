typedef struct BRDFS {

	int     total_number;
	int     used_number;
	int     *whichscene;
	int     *channel;
	int	*symmetry;
	int     normalize; /* 0-No; 1-YES */

	int	*index;
	int	*nbrdfruns;
	int	**brdfindex;
	float	**srftemp;
	float	**albedo;

	Angles  *a;

	float   ***rho_direct;
	float   ***rho_surface;
	float   ****brdfs;

	float   **hemisphere_albedo;

} Brdf;
