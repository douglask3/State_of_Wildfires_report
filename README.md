# State of Wildfire's report

Repo for analysis and plotting code used in the State of Wildfire's report:

In the follwing dirs:
* **modules** models and submodules used to generate data
*  Global ECMWF Fire Forcast (GEFF) was used to generate fwi in Figs 7, 8, 9, 16, S24, S25, S26.
* **FWI_attribution_analysis_code**
  * Scripts used to recreate the figures for the FWI attribution (1 main figure and 2 supplementary figures)
  * fwi from era5 data generated using impacts toolbox in modules/impactstoolbox-0.1-alpha.zip (weblink here: [https://github.com/douglask3/State_of_Wildfires_report/blob/main/modules/impactstoolbox-0.1-alpha.zip](https://github.com/douglask3/State_of_Wildfires_report/blob/main/modules/impactstoolbox-0.1-alpha.zip))
* **make_inputs** - code for generating input files used by ConFire. Resultant data can be found at [https://doi.org/10.5281/zenodo.11420743](https://doi.org/10.5281/zenodo.11420743).
  * isimip.py regrids isimip3a and 3b code and JULES model output into ConFire inputs. Original data can be downloaded from [https://data.ISIMIP.org/](https://data.ISIMIP.org/)
  * nrt.r regrid data obtained from cds ([https://cds.climate.copernicus.eu/](https://cds.climate.copernicus.eu/)) to provide nrt driving data for ConFire.
  * regrid_vcf.r regrids vcf vegetation cover data ready for jules-es bias correction. vcf data is downloadable using download_regrid_VCF.py.
  * jules-es bias correction is conducted using the ibicus package, found here: https://github.com/jakobwes/State-of-Wildfires---Bias-Adjustment
* **ConFire_plotting** - scripts for producing ConFire plots:
    * plot_control_stripes.r - Figure 14
    * plot_combinedControls.r - Figure's 15, 16, S13, S14
    * plot_ConFire_hist.r - Figure 18
    * plot_ConFire_ts.r - Figure 22, Table 7
    * plot-ConFire_futureFire.r - Future 23, S17-S23

