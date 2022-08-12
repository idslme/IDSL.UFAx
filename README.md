# IDSL.UFAx<img src='UFAx_educational_files/Figures/IDSL.UFAx-logo.png' width="250px" align="right" />

<!-- badges: start -->
[![Maintainer](https://img.shields.io/badge/maintainer-Sadjad_Fakouri_Baygi-blue)](https://github.com/sajfb)
[![CRAN status](https://www.r-pkg.org/badges/version/IDSL.UFAx)](https://cran.r-project.org/package=IDSL.UFAx)
![](http://cranlogs.r-pkg.org/badges/IDSL.UFAx?color=orange)
![](http://cranlogs.r-pkg.org/badges/grand-total/IDSL.UFAx?color=brightgreen)
[![Dependencies](https://tinyverse.netlify.com/badge/IDSL.UFAx)](https://cran.r-project.org/package=IDSL.UFAx)
<!-- badges: end -->

A pipeline to annotate a number of peaks from the IDSL.IPA peaklists using an exhaustive chemical enumeration-based approach. This package can perform elemental composition calculations using following 15 elements : C, B, Br, Cl, K, S, Se, Si, N, H, As, F, I, Na, O, and P.

Visit https://ufa.idsl.me/enumerating-chemical-space/exhaustive-enumeration for the detailed documentation and tutorial.

# Note
IDSL.UFAx package has a dependency on RcppAlgos and gmp R packages. In some instance to install these two packages, you may need to run the following command on a linux terminal

	sudo apt-get install libgmp-dev