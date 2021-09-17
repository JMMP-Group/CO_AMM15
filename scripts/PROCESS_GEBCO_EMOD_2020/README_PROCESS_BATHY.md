
# Here keep a record of various pre processing steps done to the Raw GEBCO/ EMODNET bathy

The scripts are pretty small for the most part typically:
   * bathy in -> process it in some way -> bathy out

There are few things that needs doing though.


   1.  Get the raw data and cut out a section for NWS, note for EMODNET some extra steps required for the tile format it comes in
   1.  We create cubes of the same data which allows us to use iris/cartopy to do the interpolation 
   1.  Because we later want to pre-smooth the bathy wit the Smagorinsky smoother this needs data that goes beyond the AMM15 domain:
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
      * Use MBells modifications of the CDFTOOLS to run the Smagorinsky Smoother
      * Requires pre formatting the raw bathy data into a format that CDF smoother expects
   1. If we want to retain  NICOS baltic we need to slice on a section for the Baltic bdy that is straight from his bathy data source
   

   
   
