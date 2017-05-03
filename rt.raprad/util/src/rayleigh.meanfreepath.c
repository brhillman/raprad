double indexrefraction_real_air();
double rayleigh_crosssection();

double rayleigh_meanfreepath(lambda, Na)
  double
    lambda,
    Na;
{
  double
    mr,
    sigma;

  mr = indexrefraction_real_air(lambda);
  sigma = rayleigh_crosssection(lambda, mr);

  return ((double) ((1./sigma*Na)));

}
