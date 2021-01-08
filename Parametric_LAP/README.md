# Parametric Multiscale LAP
 Time-varying delay estimation using local all-pass with parametric fitting

Using an all-pass filter we estimate a per sample time delay between two signals and then use a parametric fitting to give a smoother estimate of the delay

- _MultiScale_LAP_ allows the use of multiple LAP filters of different sizes in order to give a more accurate estimate then performs a parametric fitting to the delay estimate using a linear polynomial model
  - Any parametric model can be used in this case a linear polynomial basis built from Chebyshev polynomials is implemented
