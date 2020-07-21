# s_deposition
Code for the kriging analysis and modelling of atmospheric S deposition (Hinckley_et_al_NGS_2020_kriging_code.R)

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
NADP wet deposition data (http://nadp.slh.wisc.edu/ntn/) (sample csv file: NTN-All-cy.csv)
location of wet deposition sites (http://nadp.slh.wisc.edu/data/sites/NTN/?net=NTN) (sample csv file: NTNsites.csv)
CASTNET dry deposition data (https://java.epa.gov/castnet/clearsession.do) (sample csv file: Dry_Deposition_Annual_all_US_kgha_2017.csv)
location of dry deposition sites (https://www.epa.gov/castnet/castnet-site-locations) (sample csv file: castnet_locations_entire_US.csv)
PRISM precipitation quantity (ftp://prism.oregonstate.edu/monthly/ppt/) (sample files in two zipped folders: PRISM_2001_2017.zip and PRISM_1989_2000.zip)

Outputs ---------------------------------------------------------------
For each running year the code generates three GeoTiff maps including wet, dry and total S deposition for the entire US.
as example output maps for years 1989 and 2017 are included here.(in output_maps.zip)

