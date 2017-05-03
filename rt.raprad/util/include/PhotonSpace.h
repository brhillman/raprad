typedef struct Photon_Space {

	int
	  cnumber;

	double
	  zm,
	  zp;

	int
	  gridnum,
	  gridnumf;

	int
	   rt_number,
	  *rtflag,
	  *rtflagf,
	  *rtflagfwhich;

	double
	  *dz,
	  *grid,
	  *gridf;

} PhotonSpace;
