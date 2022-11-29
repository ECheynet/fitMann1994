# Fitting the Uniform Shear Model to real data
The parameters of the Uniform Shear Model [1] are estimated in the least-square sense

[![View Fitting the Uniform Shear Model to real data on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/73126-fitting-the-uniform-shear-model-to-real-data)
[![DOI](https://zenodo.org/badge/249147947.svg)](https://zenodo.org/badge/latestdoi/249147947)

## Summary
The function "fitMann" is used to estimate the parameters of the uniform shear model [1] in the least-square sense. As an example, the fitting algorithm is applied to the Great Belt dataset used by Mann [1] and the Kaimal spectral model [2,3]. The present Matlab implementation was applied for offshore wind data [4]. 

## Content

The submission contains:
- the function fitMann.m
- An interactive example file Example.mlx
- An example where fitMann is applied on incomplete data.
- a data file greatbeltData that contains digitalized spectral estimate from the Great belt experiment
- The function mannTurb.m introduced in [5]


Any comment, suggestion or question is welcomed.


## References

[1] Mann, J. (1994). The spatial structure of neutral atmospheric surface-layer turbulence. Journal of fluid mechanics, 273, 141-168.

[2] Kaimal, J. C., Wyngaard, J. C. J., Izumi, Y., & Coté, O. R. (1972). Spectral characteristics of surface‐layer turbulence. Quarterly Journal of the Royal Meteorological Society, 98(417), 563-589.

[3] Mann, J. (1998). Wind field simulation. Probabilistic engineering mechanics, 13(4), 269-282.

[4] Cheynet, E. (2019). Influence of the Measurement Height on the Vertical Coherence of Natural Wind. In Conference of the Italian Association for Wind Engineering (pp. 207-221). Springer, Cham.

[5] https://se.mathworks.com/matlabcentral/fileexchange/67055-uniform-shear-model-mann-1994
