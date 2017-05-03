typedef struct Photon_Partition {

	int
		numintervals,
		*interval,
		*runsperinterval,
		*smindex,

		*addscatterflag,
		*bbradflag,
		*bbradmaxlayers,
		*numrefl;

	float
		*wavelength,
                **bandwidth;

}  PhotonPartition;
