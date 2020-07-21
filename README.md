# s_deposition
Code for the kriging analysis and modelling of atmospheric S deposition

Description ---------------------------------------------------------
This repository includes the kriging model and input files used for E-L. S. Hinckley,
J. T. Crawford1, H. Fakhraei and C. T. Driscoll for A shift in sulfur-cycle manipulation from
atmospheric emissions to agricultural additions, Nature Geoscience 2020. This is a kriging model to
interpolate atmospheric deposition through no measured locations throughout the US using point observation
from two monitoring networks (i.e., NADP and CASTNET)
You are welcome to use or adapt this code providing you cite
E-L.S. Hinckley, J. T. Crawford1, H. Fakhraei and C.T. Driscoll for A shift in sulfur-cycle manipulation from
atmospheric emissions to agricultural additions, Nature Geoscience 2020.
To run this code a set of input data is required to be downloaded from the public domain resources.
Here, we included these datasets until year 2017 but the user will be able to download the data and
perform kriging as the datasets are getting updated.

Inputs ---------------------------------------------------------------
NADP wet deposition data (http://nadp.slh.wisc.edu/ntn/)
location of wet deposition sites (http://nadp.slh.wisc.edu/data/sites/NTN/?net=NTN)
CASTNET dry deposition data (https://java.epa.gov/castnet/clearsession.do)
location of dry deposition sites (https://www.epa.gov/castnet/castnet-site-locations)
PRISM precipitation quantity (ftp://prism.oregonstate.edu/monthly/ppt/)

Outputs ---------------------------------------------------------------
For each running year the code generates three GeoTiff maps including wet, dry and total S deposition for the entire US.
as example output maps for years 1989 and 2017 are included here.

