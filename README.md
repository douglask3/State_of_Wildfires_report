# State of Wildfire's report

Repo for analysis and plotting code used in the State of Wildfire's report:

In the follwing dirs:
* **make_inputs** - code for generating input files used by ConFire. Resultant data can be found at https://doi.org/10.5281/zenodo.11420743.
  * isimip.py regrids isimip3a and 3b code and JULES model output into ConFire inputs. Original data can be downloaded from https://data.ISIMIP.org/
  * nrt.r regrid data obtained from cds (https://cds.climate.copernicus.eu/) to provide nrt driving data for ConFire.
  * regrid_vcf.r regrids vcf vegetation cover data ready for jules-es bias correction. vcf data is downloadable using download_regrid_VCF.py.
  * jules-es bias correction is conducted using the ibicus package, found here: https://github.com/jakobwes/State-of-Wildfires---Bias-Adjustment
* **FWI_attribution_analysis_code** - Scripts used to recreate the figures for the FWI attribution (1 main figure and 2 supplementary figures)
* [[add other dirs here]]

