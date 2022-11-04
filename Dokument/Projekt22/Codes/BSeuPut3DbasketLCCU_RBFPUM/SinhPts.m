function xc = SinhPts(x_min,x_max,N,K,ell)
  xi0 = asinh((x_min-K)/ell);
  xiN = asinh((x_max-K)/ell);
  xi = linspace(xi0,xiN,N)';
  xc = K + ell*sinh(xi);