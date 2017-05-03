
double indexrefraction_real_air(lambda)
  double
    lambda;
{
  double
    mr,
    linv;

  linv = 1./lambda;

  mr  = 1. + 8342.13E-08;
  mr += (2.406030E-02 / (130. - linv*linv));
  mr += (  1.5997E-04 / (38.9 - linv*linv));

  return ((double) mr);

}
