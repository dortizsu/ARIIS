# ARIIS
This is a maintained repository for ARIIS (algorithm for robust identification of the inertial subrange).

ARIIS is an algorithm for emprically identifying Kolmogorov's inertial subrange and its spectral slope (e.g. -5/3) within
a given variance spectrum. Please see Ortiz-Suslow, Wang, Kalogiros, and Yamaguchi "A Method for Identifying Kolmogorov's
Inertial Subrange in the Velocity Variance Spectrum", Journal of Atmospheric and Oceanic Technology, 2019, for details. If you
find ARIIS helpful, please cite this paper.

This algorithim is provided as-is and open for general use, please read included License. ARIIS is actively maintained and
a component of active research. Please check repository for updated versions and please direct all queries, issues, bug 
reports to Dave Ortiz-Suslow, see contact info below. 

Contact Info: dortizsu [at] nps [dot] edu
website: dortizsuslow [dot] weebly [dot] com

Contents:
[1]ariis --> MATLAB-readable code of the ARIIS algorithm. An all-incluside-one-stop-shop program;
             feed it what it wants and it wil go.
[2]ariis_demo --> provides a simple example of how ARIIS runs and the output it gives.
[3]ariis_test.mat --> matlab archive used in demo.
[4]ariis_license

Caveats:
* ARIIS can handle temperature and concentration (e.g. water vapor) input, assuming inertial subrange in the kinetic energy has the same bandwidth in the inertial-convective and inertial-scalar subranges (Tennekes & Lumley 1972). 
* If using sonic anemometry, can handle applying path-averaging correction for CSAT-type and Solent-type sonic anemometer heads. If using open-path gas analyzers, user must input the path length and it is assumed that path is (approx.) perpindicular to the mean wind vector.
* ARIIS was designed with environmental turbulence measurements in mind (noisy turbulence measurements from an active acoustic sensor). However, the methodology is fundamentally general and can be applied to any 3D velocity and/or turbulent scalar field. Check INPUTS & CONSTANTS structures carefully.

