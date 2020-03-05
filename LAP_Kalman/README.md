# LAP + Kalman
 Time-varying delay estimation using local all-pass and Kalman filters

Using an all-pass filter we estimate a per sample time delay between two signals and then use a Kalman filter to give a smoother estimate of the delay

- _LAP_Kalman_ estimates a single (time-varying) delay between two channels using a LAP filter
   - This process can be repeated for multiple different filter sizes
   - The delay estimated from each LAP filter is assumed to be a noisy estimate of the true delay and is then passed as the input to a Kalman filter
   - As well as individual LAP+Kalman filter outputs the function also provides a fused output which takes the individual LAP outputs and performs meausrement fusion using a single Kalman filter

- _Delay_Est_Kalman_ allows comparison of the LAP, LAP+Kalman and fused LAP+Kalman using the data generated in _Signal_Generation_

## References
 APSIPA_2019 [Fast & Efficient Delay Estimation Using Local All-Pass & Kalman Filters](https://beteje.github.io/assets/pdf/2019_APSIPA.pdf)
