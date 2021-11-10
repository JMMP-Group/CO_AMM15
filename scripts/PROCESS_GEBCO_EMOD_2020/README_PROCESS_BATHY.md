
# Here keep a record of various pre processing steps done to the Raw GEBCO/ EMODNET bathy

The scripts are pretty smal for the most part typicaly:
   * bathy in -> process it in some way -> bathy out

There are few things that needs doing though.


   1.  Get the raw data and cut out a section for NWS, note for EMODNET some extra steps required for the tile format it comes in
   1.  We create cubes of the same data whcih alllows us to use iris/cartopy to do the interpolation 
   1.  Because we later want to presmooth the bathy wit the smagorinsky smoother this needs data that goes beyond the AMM15 domain:
      1. Thsu we need to generate a target grid and coordinates giles the has an extent greater than the AMM15
   1. As the data sets are referenced again LAT we need a strategy to undo that. For not that involves
      1. Using a combination of CS3X and CS20 
         * CS3X has a larger (off shelf ) coverage then CS20
         * CS3X is coarses
         * CS20 was based on an old POLCOMS run and has problems near Holland (Dykes etc not done properly)
         * Neither cover the extended AMM15, so for data beyond AMM15 we just live withouth the LAT correction
   1. We also do asplice of data for merge dataset. 
      * This is because GEBCO has less suprious features in the deep then EMODNET
      * Tidal tests seesm to show tha EMODNET does a better job of Dogger bank
      * At very shallow areas the extra detail of EMODNET actually is problematic at AMM15 1/5km resolution and we stick with the smooother GEBCO here

   
   
