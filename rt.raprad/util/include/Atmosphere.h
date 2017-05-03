typedef struct Atmosphere_Profile {
  int    numlayers;           /* Input - No of thermal layers                 */
  int    numgases;            /* Input - No of gases in atmosphere            */
  double *altitude;
  double *pressure;           /* Input - Array of pressures at each layer     */
  double *temperature;        /* Input - Array of temperatures at each layer  */
  double *airdens;            /* Input - Array of densities at each layer     */
  double *dryairdens;
  double *airdensavg;
  double *airmass;
  double *densh2o;
  double *densco2;
  double *denso3;
  double *denso2;
  double *densn2o;
  double *densch4;
  double *mxratioh2o;
  double *mxratioco2;
  double *mxratioo3;
  double *mxratioo2;
  double *mxration2o;
  double *mxratioch4;
  double *uh2oprecm;
  double *uco2amagats;
  double *uo3amagats;
  double *uo2amagats;
  double *uh2o;
  double *uco2;
  double *uo3;
  double *uo2;
  double *un2o;
  double *uch4;
} Atmosphere;
