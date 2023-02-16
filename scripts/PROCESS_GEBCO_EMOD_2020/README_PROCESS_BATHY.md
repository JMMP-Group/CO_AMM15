# Table of contents

- [Table of contents](#table-of-contents)
- [Here keep a record of various pre processing steps done to the Raw GEBCO/ EMODNET bathy](#here-keep-a-record-of-various-pre-processing-steps-done-to-the-raw-gebco-emodnet-bathy)
- [Retrieve and process source data](#retrieve-and-process-source-data)
  - [Processing GEBCO](#processing-gebco)
    - [Cut out a domain just for the North West shelf](#cut-out-a-domain-just-for-the-north-west-shelf)
    - [Convert GEBCO to Cube](#convert-gebco-to-cube)
    - [Interpolate GEBCO cube to AMM15 extended domain](#interpolate-gebco-cube-to-amm15-extended-domain)
    - [Correct for LAT and apply the operational existing AMM15 mask to the inner domain](#correct-for-lat-and-apply-the-operational-existing-amm15-mask-to-the-inner-domain)
      - [First Create a file to make the LAT correction with using Surge model data as a proxy for LAT 2D field](#first-create-a-file-to-make-the-lat-correction-with-using-surge-model-data-as-a-proxy-for-lat-2d-field)
      - [Apply the correction data to the GEBCO data on the AMM15 grid](#apply-the-correction-data-to-the-gebco-data-on-the-amm15-grid)
  - [Processing EMODNET](#processing-emodnet)
    - [Joining the Tiles](#joining-the-tiles)
    - [Removing the overlaps and creating a cube version of EMODNET data](#removing-the-overlaps-and-creating-a-cube-version-of-emodnet-data)
    - [Re-grid to extended AMM15 domain](#re-grid-to-extended-amm15-domain)
    - [Correct for LAT](#correct-for-lat)
    - [Merge EMDONET and GEBCO AMM15 data into one dataset](#merge-emdonet-and-gebco-amm15-data-into-one-dataset)
    - [Optional Pre-smooth of the data](#optional-presmooth-of-the-data)
      - [Do the smoothing](#do-the-smoothing)
  - [Cut out domain to Exact extent of AMM15 and optionally make outer 4 points of bdy equal to the outer rim](#cut-out-domain-to-exact-extent-of-amm15-and-optionally-make-outer-4-points-of-bdy-equal-to-the-outer-rim)
  - [Modify Baltic bathy](#modify-baltic-bathy)
    - [Pre Smooth the Baltic](#pre-smooth-the-baltic)
    - [Apply Nico Baltic Rim exactly (Optional)](#apply-nico-baltic-rim-exactly-optional)
  - [Merge File with Baltic modifications](#merge-file-with-baltic-modifications)

# Here keep a record of various pre processing steps done to the Raw GEBCO/ EMODNET bathy

Some Required or intermediate inputs are placed on jasmin under:

- /gws/nopw/j04/jmmp_collab/AMM15/EMODNET_GEBCO_2020/REQUIRED_INPUTS

The scripts are pretty small for the most part typically:

- bathy in -> process it in some way -> bathy out

There are few things that needs doing though.

   1. Get the raw data and cut out a section for NWS, note for EMODNET some extra steps required for the tile format it comes in
   1. We create cubes of the same data which allows us to use iris/cartopy to do the interpolation
   1. Because we later want to pre-smooth the bathy with the Shapiro smoother this needs data that goes beyond the AMM15 domain:
      1. Thus we need to generate a target grid and coordinates files the has an extent greater than the AMM15
   1. As the data sets are referenced against LAT we need a strategy to undo that. For now that involves
      1. Using a combination of CS3X and CS20
         - CS3X has a larger (off shelf ) coverage then CS20
         - CS3X is coarse
         - CS20 was based on an old POLCOMS run and has problems near Holland (Dykes etc not done properly)
         - Neither cover the extended AMM15, so for data beyond AMM15 we just live without the LAT correction
   1. We also do a splice of data for merged dataset.
      - This is because GEBCO has less spurious features in the deep then EMODNET
      - Tidal tests seems to show that EMODNET does a better job of Dogger bank
      - At very shallow areas the extra detail of EMODNET actually is problematic at AMM15 1.5 km resolution and we stick with the smoother GEBCO here

   1. Pre-Smooth use CDFTOOLS
      - We need the larger than AMM15 domain as input
      - Use MBells modifications of the CDFTOOLS to run the Shapiro Smoother
      - Requires pre formatting the raw bathy data into a format that CDF smoother expects
   1. If we want to retain NICOS Baltic we need to slice on a section for the Baltic bdy that is straight from his bathy data source

# Retrieve and process source data

The source data considered is that of EMODNET and GEBCO 2020.

- GEBCO:
  - <https://www.gebco.net/data_and_products/gridded_bathymetry_data/#global>  (this at time of writing is  now up to 2021)
- EMODNET:
  - <https://portal.emodnet-bathymetry.eu/>   (This comes in a selection of tiles with some overlap at the edges that needs care of when processing)

## Processing GEBCO

### Cut out a domain just for the North West shelf

Sample files on JASMIN:

- /gws/nopw/j04/jmmp_collab/AMM15/EMODNET_GEBCO_2020/GEBCO_2020/SOURCE

They are the global data set of the elevations and also contain the LSM in the tid file.

We make a more manageable cut out of each for the NWS:

```bash
 ncks -d lon,33000,52000 -d lat,29000,40000 GEBCO_2020.nc NWS_CUT_GEBCO.nc
 ncks -d lon,33000,52000 -d lat,29000,40000 GEBCO_2020_TID.nc NWS_CUT_GEBCO_2020_TID.nc
```

With a copy on JASMIN here:

- /gws/nopw/j04/jmmp_collab/AMM15/EMODNET_GEBCO_2020/GEBCO_2020

The data needs to be mapped onto the AMM15 grid. To do so we can make use of the iris interpolation.

### Convert GEBCO to Cube

In order to use the interpolator we convert the GEBCO data into an iris cube. To do that use:

- [GEBCO_PROCESS/MAKE_GEBCO_CUBE.py](GEBCO_PROCESS/MAKE_GEBCO_CUBE.py)

This takes as input the path to the cutout of the NWS GEBCO data (NWS_CUT_GEBCO.nc) and the path
of the dir to store the resultant cube of data (GEBCO_CUBE.nc).

### Interpolate GEBCO cube to AMM15 extended domain

In the later stages when we use the Shapiro smoother we need data that goes beyond AMM15 boundaries.
Thus we actually map the GEBCO data onto a cube who's lat lon extents are beyond the AMM15 lat lon extents,
but use the same underlying grid.

In order to expand the grid we need to obtain the underlying AMM15 grid to begin with.
This is defined in the grid file AMM15_ROTATED_CS.nc and stored on JASMIN under:

- wget <https://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/EMODNET_GEBCO_2020/REQUIRED_INPUTS/AMM15_ROTATED_CS.nc>
- wget <https://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/EMODNET_GEBCO_2020/GEBCO_2020/NWS_CUT_GEBCO_2020_TID.nc>

We use the script

- GEBCO_PROCESS/EXPAND_AMM15_CUBE.py

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

This script also computes a NEMO style coordinates file for the extended AMM15 routine using the function.

- output_nemo_coords
From
- EXPAND_AMM15_CUBE.py

```
the following arguments are required: -a/--AMM15_PATH, -i/--INLSM_DIR, -c/--INCUBE_DIR, -o/--OUT_DIR -y/--YEAR
```

e.g.:

```
python EXPAND_AMM15_CUBE.py -a ~/PATH_TO_COORDS/AMM15_ROTATED_CS.nc -i PATH_TO_GEBCO_LSM -c PATH_TO_GEBCO_CUBE -o PATH_TO_OUTPUT -y GEBCO_VERSION_YEAR
```

We process the data outside the inner core AMM15 domain and inside differently.

- We compute a version of the bathymetry that is extrapolated everywhere regardless of LSM.
- We compute a version of the bathymetry that is masked by the GEBCO mask.
- We merge the two but retain the extrapolated version in the inner domain and the LSM in the outer domain.
  - Later we apply the operational LSM on the inner domain ahead of Shapiro smoothing when we also apply the LAT.

### Correct for LAT and apply the operational existing AMM15 mask to the inner domain

EMODNET is explicitly quoted as being referenced to Lowest Astronomical Tide (LAT).
When we compare EMODNET and GEBCO data it appears that there are regions where they are basically the same
data as EMODNET defaults to GEBCO as its basemap. Thus it is reasonable to assume that GEBCO  must be LAT referenced.
Thus we need to apply an LAT correction to this data consistently in the way we intend to with EMODNET.

#### First Create a file to make the LAT correction with using Surge model data as a proxy for LAT 2D field

In order to do that we need an approximation of the LAT. Our current estimation relies on  surge models CS3X and CS20.

CS3X covers the AMM15 domain and CS20 is at a higher resolution but only covers out to the shelf break.

Our approach is to merge the two solutions and to use the CS3X solution beyond the shelf break.

The data for CS3X and CS20 is kindly provided by Colin Bell. It comes in ASCII format and in converted
into mask netcdf grids using the routines:

- convertCS3X.py
  - converts:  CS3X_stat.txt -> CS3X_stat.nc
- convertCS20.py
  - converts:   cs20_stat.txt  -> cs20_stat.nc

Note Jenny Graham already processed the CS3X data so we can use that file directly instead.

Once we have both netcdf files as can read them in and combine them together on the AMM15 grid.
This is done with:

- SURGE_LAT_CORRECTION/Combine_Surge.py

It takes as an input argument -i the path for the location of the required surge netcdf files and the AMM15 grid file,
and an argument -o where to output the data to.

Note can retrieve the merged file :

`
wget https://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/EMODNET_GEBCO_2020/REQUIRED_INPUTS/MERGE_JENNY_C3X_COLIN_CS20.nc
`

#### Apply the correction data to the GEBCO data on the AMM15 grid

With the merge CS3X CS20 LAT proxy field we can correct the inner domain of the AMM15 (for where we have LAT data) using the script:

- GEBCO_PROCESS/Correct_LAT_apply_op_mask.py

This script takes a number of inputs:

- OP_LSM :
  - The location of the operational LSM e.g. EMODNET_LSM_v2.nc (created by J graham from EMODNET to AMM15 plus fill in lakes etc)
    - wget <https://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/EMODNET_GEBCO_2020/REQUIRED_INPUTS/EMODNET_LSM_v2.nc>
- CS3X_CS20:
  - The location of the merged surge data for LAT correction (valid only on AMM15 inner domain)
    - wget <https://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/EMODNET_GEBCO_2020/REQUIRED_INPUTS/MERGE_JENNY_C3X_COLIN_CS20.nc>
- BATHY_DATA:
  - The location the GEBCO data on the extended AMM15 grid

    - LAT_LON:
      - The location of the file containing the extended AMM15 domain lat lon grid

- OUT_FILE:
  - Where to write the final output file after processing

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

- bash runmerge.sh -i  InputDirPath -o OutputDirPath -p yourPythonCommand

The input path should contain the EMODNET tiles and then it loops over C to F calling *merge_xarray.py* to create rows of data.
Then it does a final merge of the rows into a unified block of data using *final_merge.py*

The resultant merged file is ALLmerge.nc

### Removing the overlaps and creating a cube version of EMODNET data

The above data is stripped of overlaps along the edges and then put into cube format by

- *MAKE_EMODNET_CUBE.py*

Typical usage is:

- python3.8 MAKE_EMODNET_CUBE.py -i PathToALLmerge.nc/ -o OutPutDir/

To reduce RAM usage this uses xarray chunking.

```python
EMODNET_RAWB = xr.open_mfdataset( '{}/ALLmerge.nc'.format(args.IN_DIR[0]), parallel=True,chunks=({"lat": 100, "lon": -1 }) )
```

It searches the lats and lons for overlaps and creates an array of indexes so that the global array can be indexed without the overlaps,
Note in some cases the indexes are almost the same but only differ by a tiny amount so we exclude those cases too.

```python
ind_x = np.zeros(1,dtype=int)
i=1
while(i < np.size(EMODNET_LON.data[:]) ):
    if((EMODNET_LON.data[i] - EMODNET_LON.data[i-1])<1e-10):
        k = i-1
        print( "caught a repeat", EMODNET_LON.data[i],EMODNET_LON.data[i-1])
        while( (EMODNET_LON.data[i] - EMODNET_LON.data[k]) < 1e-10):
            print( "Caught a further repeat", EMODNET_LON[i].data,EMODNET_LON.data[k],-EMODNET_LON.data[i] , EMODNET_LON.data[k])
            i=i+1
    ind_x = np.append(ind_x[:],i)
    i=i+1
print(ind_x)

ind_x = xr.DataArray(ind_x, dims=["lon"])
```

similar for y and then index the original data_and_products

```python
dask.config.set(**{'array.slicing.split_large_chunks': False})
EMODNET_BATHY = EMODNET_BATHY[ind_y,ind_x]
```

We rename the data variable for consistency with later scripts

```python
EMODNET_BATHY = EMODNET_BATHY.rename("sea_floor_depth_below_geoid")
```

And convert to a cube

```python
EMODNET_cube = EMODNET_BATHY.to_iris()
```

We also want the same coordinate system as AMM15 for re-gridding later.

The output

```python
iris.save(EMODNET_cube, '{}/EMODNET_v2020_NO_REPEAT_LAT_LON.nc'.format(args.OUT_DIR[0]))
```

is used in the next stage where we re-grid onto the extended AMM15 domain

### Re-grid to extended AMM15 domain

To regrid to the extended AMM15 domain we use the script:

- EXPAND_AMM15_CUBE.py
This is amended from the similar version as in GEBCO as we ran into RAM limitations.

The iris re-grid function seems to grab all the data despite using lazy loading etc.
So to over come this we break the target grid into small sections. We also choose
src grid that covers the target grid and an extra wind of say 500 points only.
This reduces the RAM required. The user can split the domain up into smaller sections
if they have limited RAM. This is changed by setting sub_size(currently hard-coded)

```python
sub_size  = 300
```

Arguments passed are:

- -a ~/EMODNET_GEBCO_2020/REQUIRED_INPUTS/AMM15_ROTATED_CS.nc
  - That is the path to the example AMM15 rotated grid
- -c ~/EMODNET_GEBCO_2020/REQUIRED_INPUTS  
  - the path to where the EMODNET src data has been mapped to a cube
  - e.g. EMODNET_v2020_NO_REPEAT_LAT_LON.nc
- -o  ~/EMODNET_GEBCO_2020/REQUIRED_INPUTS
  - The path where we want to store our oututs:
    - MASK_EMODNET_vDec2020_ON_EXPAND_AMM15.nc
    - EXTRAPOLATE_EMODNET_vDec2020_ON_EXPAND_AMM15.nc
    - **MASK_EXTRAPOLATE_EMODNET_vDec2020_ON_EXPAND_AMM15.nc**

When the code runs it calculates the number of blocks or sections required based on sub_size.
It creates a SUBSECTION subdir into which is writes each subsection.
To get the src data correct it finds the  lat lon box that encompasses the target
section and adds on extra perimeter so that filling in from nearest valid sea point is consistent.
e.g.

```python
find closest indices
# Min
  value,lat_idx = find_nearest(EMODNET_RAW_cube.coord('latitude').points[:], np.min(lats ))
  lat_min_idx = lat_idx - 500 # ( for safety for fill interp)

```

etc. and then the subsection of src data is:

```python
  subsection_emodnet_raw_cube = EMODNET_RAW_cube[lat_min_idx:lat_max_idx, lon_min_idx:lon_max_idx ]
```

As in GEBCO case it nans fills masks etc as before.
To create the final joined up data set it reads in all the subsections using open_mfdataset:

```python
EXTRAPOLATE_EMODNET_ON_EXPANDAMM15 = xr.open_mfdataset(
        '{}/SUBSECTION/00???_FILL_SUBSECTION_CUBE.nc'.format(args.OUT_DIR[0]),
        combine = 'nested',
        concat_dim = 'grid_latitude',
        parallel=True)

```

and then write back out as a unified dataset e.g.:

```python
with ProgressBar():
  MASK_EXTRAPOLATE.to_netcdf('{}/MASK_EXTRAPOLATE_EMODNET_vDec2020_ON_EXPAND_AMM15.nc'.format(args.OUT_DIR[0]))
```

### Correct for LAT

We ended up using the merge of Jenny's CS3X data not Colin's CS3X Data
So the above section needs more information

The script we use to apply the LAT corrections and the operational mask is
`Correct_LAT_Apply_op_mask.py`

```bash
python Correct_LAT_apply_op_mask.py -o .../EMODNET_GEBCO_2020/REQUIRED_INPUTS/EMODNET_LSM_v2.nc -c .../EMODNET_GEBCO_2020/REQUIRED_INPUTS/MERGE_CS3X_COLIN_CS20.nc  -b  ..../GEB_EMODNET_PROCESSING_GIT_VERSION/MASK_EXTRAPOLATE_EMODNET_vDec2020_ON_EXPAND_AMM15.nc .../GEB_EMODNET_PROCESSING_GIT_VERSION/expand_amm15.coordinates.nc -f .../GEB_EMODNET_PROCESSING_GIT_VERSION/FINAL_EMODNET_LAT_CORRECTED_EXPANDED_AMM15.nc
```

### Merge EMDONET and GEBCO AMM15 data into one dataset

Once we have EMODNET and GEBCO with Correction to LAT on the extended AMM15 grid we can merge the two.

(see further notes on why in README_MERGE_EXPANDED.md)

Basic idea is simple enough with GEBCO in the deep, EMODNET 100-5 m and GEBCO to the coast

We use

```bash
python3.8 CORRECTED_EMDONET_INPUT_EXPANDED_MERGE_GEBCO_DEEP_to200-100M_EMODNET_TO_10-5M_GEBCO_TO_COAST.py -e ~/scratch/FINAL_EMODNET_LAT_CORRECTED_EXPANDED_AMM15.nc -g ~/scratch/FINAL_GEBCO_LAT_CORRECTED_EXPANDED_AMM15.nc  -c ~/scratch/expand_amm15.coordinates.nc -o ~/scratch/
```

where

- -e/--EMOD_FILE :
  - The location of the EMODNET data with LAT correction on extended AMM15 grid
- -g/--GEB_FILE :
  - The location of the GEBCO data with LAT correction on extended AMM15 grid
- -c/--COORD_FILE :
  - The location of the Extended AMM15 coordinates file
- -o/--OUT_DIR :
  - The path to the output directory

### Optional Presmooth of the data

N.B.: this PRESMOOTH step is not actually used in Enda's bathymetry since sensitivity test showed that tides were degraded. 

We can use Shapiro pre smoothing of the data. We expanded the domain to allow this to happen,
as at the boundary edge the smoother will not have data to work from. By making it larger
than AMM15 means this edge artefact will remain outside the final AMM15 domain.

see Branch  Mikes_Shapiro_Changes of git@github.com:endaodea/CDFTOOLS_4.0_ISF.git

```bash
git clone --branch Mikes_Shapiro_Changes git@github.com:endaodea/CDFTOOLS_4.0_ISF.git
```

Note to compile it needed all of netcdf

```bash
sudo  apt-get install *netcdf*
```

create make.macro in src :

```bash
# Makefile for CDFTOOLS
#    $Rev: 522 $
#    $Date: 2011-06-17 12:50:13 +0200 (Fri, 17 Jun 2011) $
# --------------------------------------------------------------
#
#NETCDF_ROOT=/opt/netcdf-4.1.1-gfortran
NETCDF_ROOT=/usr
NCDF = -I$(NETCDF_ROOT)/include -L$(NETCDF_ROOT)/lib -lnetcdff -lnetcdf



#F90=gfortran -v
F90=gfortran
MPF90=gfortran
#OMP=-fopenmp
NC4 = -D key_netcdf4 -D key_CMIP6
FFLAGS= -O0  $(NCDF) $(NC4) $(CMIP6)   -ffree-line-length-512  -fconvert=big-endian
LMPI=-lmpich

INSTALL=$(HOME)/bin
INSTALL_MAN=$(HOME)/man
OMP=-openmp
```

then make

Need to pre-format the data, use:

```bash
python xarray_format_notemplate.py
```

Note change the input data as needed, results in  **TESTT.nc**

#### Do the smoothing

```bash
./bin/cdfsmooth -f TESTT.nc -c 2 -t S -npass 3 -lap T -lsm TT -nc4 
```

This will result in pre-smoothed file under:

```bash
TESTT.ncS23TT
```

Note this needs documentation final processing:

- Cut it down to the actual AMM15 domain
- Option to enforce the bathymetry is match at the outer perimeter for x points
- Dealing with the Baltic
  - Need to use Rmax = 0.1 in the Baltic region
  - Need to exactly match the Baltic bdy itself for Nicos boundary

Turns out the best bathy actually does not include the Shapiro filter but we include the above in case it is useful in the future with different vertical coordinate systems

## Cut out domain to Exact extent of AMM15 and optionally make outer 4 points of bdy equal to the outer rim

The idea of extending the domain beyond the AMM15 domain extents was to allow for the Shapiro Smoother. This smoother actual makes for worse tides but the methodology i s archived here in case of future need.

To go back from an extended grid to the actual AMM15 grid for
the merge G-E-G data we can use the script:

```
Cut_and_copy_bdy_perimeter.py
```

which requires the input file to be modified e.g.

```
CORRECTED_EXPANDED_MERGE_GEBCO_DEEP_TO_200-100_EMODNET_TO_10-5_GEBCO_TO_COAST_amm15.bathydepth.co7.cs3x.cs20.nc
```

and the output dir to put the outputted files

```
python Cut_and_copy_bdy_perimeter.py -i path_to_input/CORRECTED_EXPANDED_MERGE_GEBCO_DEEP_TO_200-100_EMODNET_TO_10-5_GEBCO_TO_COAST_amm15.bathydepth.co7.cs3x.cs20.nc -o  path_to_output
```

It outputs 2 files:

1. Just the simple cut out of the AMM15 domain
2. The same but wit the bdy perimeter modified

for the BDY COPY what it does is:

```python
for i in range(4):
    dscut[   i  ,   :  ] = dscut[   4  , :  ]
    dscut[   i  ,   :  ] = dscut[   4  , :  ]
    dscut[   :  ,   i  ] = dscut[   :  , 4  ]
    dscut[ -1-i ,   :  ] = dscut[   -5 , :  ]
    dscut[   :  , -1-i ] = dscut[   :  , -5 ]
```

this emulates what is done in CO7.

## Modify Baltic bathy

There are 2 aspects to this

1. Pre Smooth the bath in the Kattegat-Baltic to Rmax 0.1
2. Batch the BDY RIM to Nicos Baltic

The first part is done with pre-processing smoother

The second part currently is done by directly matching existing data.

Would be better to replace this by code.

### Pre Smooth the Baltic

This can be done with : RMAX_BATHY_LIMTED_AREA.py

```bash
python -i RMAX_BATHY_LIMTED_AREA.py -i path_to_Input_file -o path to output file 
```

the resultant file stores the original what area it applies the smoother to and an anomaly of the smoothed data compared to the original.

### Apply Nico Baltic Rim exactly (Optional)

To replicate Nicos BDY we can apply it directly to the bathymetry. This is just a simple cut and paste.

This depends on Nicos original smoothed file all we want from this is the bdy data on the rim:

- wget <https://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/EMODNET_GEBCO_2020/REQUIRED_INPUTS/NICO_MODIFIED_BATHY.nc>

Typical usage:

```bash
python splice_nic_bal_rim.py  -i path_to_/SMOOTH_BDY_COPY_CUTAMM15_CORRECTED_EXPANDED_MERGE_GEBCO_DEEP_TO_200-100_EMODNET_TO_10-5_GEBCO_TO_COAST_amm15.bathydepth.co7.cs3x.cs20.nc -n path_to_/REQUIRED_INPUTS/NICO_MODIFIED_BATHY.nc -o path_to_output 
```

## Merge File with Baltic modifications

resultant file (this has no pre smoothing done except in the Baltic region)

```bash
- wget https://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/EMODNET_GEBCO_2020/ADD_NICO_BALTIC_BDY_RIM_SMOOTH_BDY_COPY_CUTAMM15_CORRECTED_EXPANDED_MERGE_GEBCO_DEEP_TO_200-100_EMODNET_TO_10-5_GEBCO_TO_COAST_amm15.bathydepth.co7.cs3x.cs20.nc
```
