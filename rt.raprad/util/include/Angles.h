typedef struct SUN {

	int	sunz_number;

	float   sunz_theta;
	float   sunz_mu;

	int	suna_number;

	float   suna_theta;
	float   suna_mu;

	double  sun2earth_distance;

        float   gmt_time;

} Sun;

typedef struct ANGLES {

	int	 sunz_number;
	float   *sunz_theta;
	float   *sunz_mu;

	int	 suna_number;
	float   *suna_theta;
	float   *suna_mu;

	int	 z_number;
	float   *z_theta;
	float   *z_mu;
	float   *z_weights;

	int	 a_number;
	float   *a_theta;
	float	 a_dtheta;
	float   *a_mu;
	float   *a_weights;

	int	 zxa_number;

} Angles;
