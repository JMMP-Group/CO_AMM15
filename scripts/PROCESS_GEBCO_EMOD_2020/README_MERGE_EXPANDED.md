# Merging the expanded GEBCO and EMODNET data on AMM15 ahead of Smagornisky Smoothing.

### Use Atom
Ctrl+Shift+M

## EXPANDED_MERGE_GEBCO_DEEP_to200-100M_EMODNET_TO_10-5M_GEBCO_TO_COAST.py

Basic idea is to read in GEBCO and EMODNET already interpolated to the expanded AMM15 DOMAIN and then do linear interpolation


*  Read In
```python
gebfname =('../../PROCESS_GEBCO/GEBCO_2020_expand_amm15.bathydepth.co7.cs3x.cs20.nc')
emodfname =('../EXPANDED_AMM15_MERGE_GEBCO_EMODNET_10_25_2020_amm15.bathydepth.co7.cs3x.cs20.nc')
```

* Set top and bottom levels
```python
top_1 = 0200.
bot_1 = 0100.
top_2 = 0010.
bot_2 = 0005.
```

* set up linear interpolation for deep to mid waters
* This is goin to be GEBCO from deep to upper top (200 m)  then linear interpolate to EMODNET by upper bottom (100 m)


```python
ramp_pass_1    = np.zeros( np.shape(emodnet_bathy) )
ramp_pass_1[:] = 1 - (emodnet_bathy[:]-bot_1 )/top_1


ramp_pass_1[ np.where( emodnet_bathy[:] >=  top_1 )] = 0.
ramp_pass_1[ np.where( emodnet_bathy[:]  <  bot_1 )] = 1.

emodnet_bath_merge_1 =  emodnet_bathy *ramp_pass_1[:] +  (1-ramp_pass_1[:])*gebco_bathy[:]
```
Then do the next interpolation from top_2 to bot_2 (emodnet back to gebco)

```python
#run second pass for shallow merge
ramp_pass_2    = np.zeros(np.shape(emodnet_bath_merge_1))
ramp_pass_2[:] = 1 - (emodnet_bath_merge_1[:]-bot_2 )/top_2

ramp_pass_2[ np.where( emodnet_bath_merge_1[:] >=  top_2 )] = 0.
ramp_pass_2[ np.where( emodnet_bath_merge_1[:]  <  bot_2 )] = 1.


emodnet_bath_merge_2 =  gebco_bathy *ramp_pass_2[:] +  (1-ramp_pass_2[:])*emodnet_bath_merge_1[:]
```

Read in the expanded coordinates file and then write all out to file

```python
ncfile = Dataset('../EXPANDED_MERGE_GEBCO_DEEP_TO_%d-%d_EMODNET_TO_%d-%d_GEBCO_TO_COAST_amm15.bathydepth.co7.cs3x.cs20.nc'%(int(top_1),int(bot_1),int(top_2),int(bot_2)),'w')
y_dim = ncfile.createDimension('x', x)     # Y
x_dim = ncfile.createDimension('y', y)     # X
latout = ncfile.createVariable('lat', 'f4', ('y','x',))
lonout = ncfile.createVariable('lon', 'f4', ('y','x',))
bathy = ncfile.createVariable('Bathymetry', 'f4', ('y','x',))

latout[:]=lat
lonout[:]=lon

bathy[:]=emodnet_bath_merge_2
```

After that we can then look to using the Smagorinsky interpolation
