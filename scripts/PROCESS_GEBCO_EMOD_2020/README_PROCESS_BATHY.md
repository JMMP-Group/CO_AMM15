
# Here keep a record of various pre processing steps done to the Raw GEBCO/ EMODNET bathy


Some Required or intermediate inputs are placed on jasmin under:
   * /gws/nopw/j04/jmmp_collab/AMM15/EMODNET_GEBCO_2020/REQUIRED_INPUTS

The scripts are pretty small for the most part typically:
   * bathy in -> process it in some way -> bathy out

There are few things that needs doing though.


   1.  Get the raw data and cut out a section for NWS, note for EMODNET some extra steps required for the tile format it comes in
   1.  We create cubes of the same data which allows us to use iris/cartopy to do the interpolation
   1.  Because we later want to pre-smooth the bathy wit the Shapiro smoother this needs data that goes beyond the AMM15 domain:
      1. Thus we need to generate a target grid and coordinates files the has an extent greater than the AMM15
   1. As the data sets are referenced again LAT we need a strategy to undo that. For not that involves
      1. Using a combination of CS3X and CS20
         * CS3X has a larger (off shelf ) coverage then CS20
         * CS3X is coarse
         * CS20 was based on an old POLCOMS run and has problems near Holland (Dykes etc not done properly)
         * Neither cover the extended AMM15, so for data beyond AMM15 we just live without the LAT correction
   1. We also do a splice of data for merge dataset.
      * This is because GEBCO has less spurious features in the deep then EMODNET
      * Tidal tests seems to show that EMODNET does a better job of Dogger bank
      * At very shallow areas the extra detail of EMODNET actually is problematic at AMM15 1.5 km resolution and we stick with the smoother GEBCO here

   1. Pre-Smooth use CDFTOOLS
      * We need the larger then AMM15 domain as input
      * Use MBells modifications of the CDFTOOLS to run the Shapiro Smoother
      * Requires pre formatting the raw bathy data into a format that CDF smoother expects
   1. If we want to retain  NICOS Baltic we need to slice on a section for the Baltic bdy that is straight from his bathy data source




# Retrieve and process source data

The source data considered is that of EMODNET and GEBCO 2020.

  * GEBCO:
      * https://www.gebco.net/data_and_products/gridded_bathymetry_data/#global  (this at time of writing is  now up to 2021)
  * EMODNET:
      * https://portal.emodnet-bathymetry.eu/   (This comes in a selection of tiles with some overlap at the edges that needs care of when processing)

## Processing GEBCO

### Cut out a domain just for the North West shelf

Sample files on jamsin:

   *  /gws/nopw/j04/jmmp_collab/AMM15/EMODNET_GEBCO_2020/GEBCO_2020/SOURCE


They are the global data set of the elevations and also contain the LSM in the tid file.

we make a more manageable cut out of each for the NWS:

```bash
 ncks -d lon,33000,52000 -d lat,29000,40000 GEBCO_2020.nc NWS_CUT_GEBCO.nc
 ncks -d lon,33000,52000 -d lat,29000,40000 GEBCO_2020_TID.nc NWS_CUT_GEBCO_2020_TID.nc
```

With a copy on jasmin here:

   * /gws/nopw/j04/jmmp_collab/AMM15/EMODNET_GEBCO_2020/GEBCO_2020


The data needs to be mapped onto the AMM15 grid. To do we we can make use of the iris interpolation.

### Convert GEBCO to Cube

In order to use the interpolator we convert the GEBCO data into an iris cube. To do that use:

   * GEBCO_PROCESS/MAKE_GEBCO_CUBE.py

This takes as input the path to the cutout of the NWS GEBCO data (NWS_CUT_GEBCO.nc) and the path
of the dir to store the resultant cube of data (GEBCO_CUBE.nc).


### Interpolate GEBCO cube to AMM15 extended domain

In the later stages when we use the Shapiro smoother we need data that goes beyond AMM15 boundaries.
Thus we actually map the GEBCO data onto a cube whos lat lon extents are beyond the AMM15 lat lon extents,
but use the same underlying grid.

In order to expand the grid we need to obtain the underlying AMM15 grid to begin with.
This is defined in the grid file AMM15_ROTATED_CS.nc and stored on jasmin under:

   * /gws/nopw/j04/jmmp_collab/AMM15/EMODNET_GEBCO_2020/AMM15_ROTATED_CS.nc

We use the script

   * GEBCO_PROCESS/EXPAND_AMM15_CUBE.py

to interpolate the GEBCO data onto the expanded AMM15 grid.

It takes as input the path to the AMM15 rotated CS file, the path to the GEBCO cube
and the GEBCO mask file (TID).  It also takes as an argument the path to where to store the resultant  
interpolated files.


The mask on AMM15 will be defined to be exactly that as in the operational case by reading in a reference LSM.
This will be applied when we later apply the LAT correction ahead of smoothing.

For the domain that goes beyond the AMM15, we derive the mask based on the GEBCO mask. This will have many lakes
and fine scale estuaries etc. that do not connect with the sea due to lack of required resolution at 1.5 km etc.
We don\'t at this stage want to redo the cleaning up of the mask for the inner AMM15 hence we just apply what we
already have. But we don\'t need to worry about unconnected lakes etc. in the domain beyond AMM15 as it is only used
for the Shapiro filter to have data that goes beyond the AMM15 boundaries.

This script also computes a NEMO style coordinates file for the extended AMM15 routine using the function
   * output_nemo_coords
From
   * EXPAND_AMM15_CUBE.py


We process the data outside the inner core AMM15 domain and inside differently.
   * We compute a version of the bathymetry that is extrapolated everywhere regardless of LSM
   * We compute a version of the bathymetry that is  masked by the GEBCO mask
   * we merge the two but retain the extrapolated version in the inner domain and the LSM in the outer domain
      * later we apply the operational LSM on the inner domain ahead of Shapiro smoothing when we also apply the LAT



### Correct for LAT and apply the operational existing AMM15 mask to the inner domain.


EMODNET is explicitly quoted as being references to Lowest Astronomical Tide (LAT).
When we compare EMODNET and GEBCO data it appears that there are regions where they are basically the same
data as EMODNET defaults to GEBCO as its basemap. Thus it is reasonable to assume that GEBCO to must be LAT referenced.
Thus we need to apply an LAT correction to this data consistently in the way we intend to with EMODNET.

#### First Create a file to make the LAT correction with usgin Surge model data as a proxy for LAT 2D field
In order to do that we need an approximation of the LAT. Our current estimation relies on  surge models CS3X and CS20.

CS3X covers the AMM15 domain and CS20 is at a higher resolution but only covers out to the shelf break.

Our approach is to merge the two solutions and to use the CS3X solution beyond the shelf break.

The data for CS3X and CS20 is kindly provided by Colin Bell. It comes in ASCII format and in converted
into mask netcdf grids using the routines:

   * convertCS3X.py
      * converts:  CS3X_stat.txt -> CS3X_stat.nc
   * convertCS20.py
      * converts:   cs20_stat.txt  -> cs20_stat.nc

Note Jenny Graham already processed the CS3X dat so we can use that file directly instead.



Once we have both netcdf files as can read them in and combine them together on the AMM15 grid.
This is done with:

   * SURGE_LAT_CORRECTION/Combine_Surge.py


It takes as an input argument -i the path for the location of the required surge netcdf files and the AMM15 grid file,
and an argument -o where to output the data to.




#### Apply the correction data to the GEBCO data on the AMM15 grid
With the mrege CS3X CS20 LAT proxy field we can correct the inner domain of the AMM15 (for where we have LAT data) using the script:

   * GEBCO_PROCESS/Correct_LAT_apply_op_mask.py

This script takes a number of inputs:
   * OP_LSM :
       * The location of the operational LSM e.g. EMODNET_LSM_v2.nc (create by J graham from EMODNET to AMM15 plus fill in lakes etc)
   * CS3X_CS20:
       * The location of the merged Surge data for LAT correction (valid only on AMM15 inner domain)
   * BATHY_DATA:
       * The location the GEBCO data on the extended AMM15 grid
    * LAT_LON:
       * The location odf the file containinf the extended AMM15 domain lat lon grid
   * OUT_FILE:
       * Where to write the final output file after processing
It applies the LAT correction only to the inner domain:
```python
inflate_lat = 100
inflate_lon = 100
# Thus we effectively have a core part of the domain and an outer part
input_bathy_amm15core = input_bathy[inflate_lat:-inflate_lat,inflate_lon:-inflate_lon
```



## Processing EMODNET


### Joining the Tiles
EMODNET data comes in a set of tiles. Each Tile has has an overlap with its neighbour. We need to splice all the tiles together into one dataset.
And then remove any overlaps that occurring the single unified data set.

A script is called to join the tiles typical format is :
   * bash runmerge.sh -i  InputDirPath -o OutputDirPath -p yourPythonCommand

The input path should contain the EMODNET tiles and then it loops over C to F calling *merge_xarray.py* to create rows of data.
Then it does a final merge of the rows into a unified block of data using *final_merge.py*

The resultant merged file is ALLmerge.nc

### Removing the overlaps and creating a cube version of EMODNET data

The above data is stripped of overlaps along the edges and then put into tcube format by
   * *MAKE_EMODNET_CUBE.py*

Typical usage is:
   * python3.8 MAKE_EMODNET_CUBE.py -i PathToALLmerge.nc/ -o OutPutDir/

This is quite RAM intensive.
