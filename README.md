# LAP_DelayEstimation
 Time-varying delay estimation using local all-pass filters

Using an all-pass filter we estimate the time delay between two or more signals. The local all-pass (LAP) filter framework allows estimation of time-varying delays by using a short window to estimate a per sample delay.

The code here is based on the 2D image registration code here: [https://sites.google.com/site/cwsgilliam/LAP](https://sites.google.com/site/cwsgilliam/LAP)

- _LAP_1D_ estimates a single (time-varying) delay between two or more channels using a LAP filter
- _MultiScale_LAP_ allows the use of multiple LAP filters of different sizes in order to give a more accurate estimate

  Note that for each scale this function estimates the time-varying delay for the whole signal, Gaussian smooths the estimate then aligns the signals using _imshift_ and repeats for the next scale.

- _Delay_Est_ allows comparison of the two using the data generated in _Simple_EMG_Model_.
   - The _Simple_EMG_Model_ generates multiple channels of data with a choice of different conduction velocities (and hence delays).

---
## References
 [Time-Varying Delay Estimation Using Common Local All-Pass Filters with Application to Surface Electromyography](https://beteje.github.io/assets/pdf/2018_ICASSP.pdf), ICASSP 2018