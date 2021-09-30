""" Script to produce coordinate file for amm15 setup.

 Variables needed in nemo coord file are:
 lambda (lon), phi(lat), e1 and e2 (spacing)
 for t,u,v and f points.
 See p54 of nemo book.

 Import libs for netcdf reading etc

 inputs are lat,lon and path to ourput coords file to
 outputs coordinate file

"""

def output_nemo_coords(rlat,rlon,outdir):
  import numpy as np
  import iris
  import iris.analysis.cartography

  from netCDF4 import Dataset

  from pathlib import Path
  import argparse, sys

  import sys, platform

  from datetime import datetime
  import subprocess



  nlat = np.size(rlat)
  nlon = np.size(rlon)

  print( "rlat is :",rlat)
 
  shape2d= np.shape(np.zeros((nlat,nlon)))



# NB. rlat/rlon are rotated - needed to calculate e1/e2, 
# but need unrot coord for glam/gphi.

  radius = 6371229.0  # radius of earth in meters
  radian = np.pi / 180.0

  gradlat = np.gradient(rlat)
  gradlon = np.gradient(rlon)

  lat_2d = np.ones(shape2d)
  gradlat_2d =np.ones(shape2d)

  for i in range(0,np.size(rlon)):
    lat_2d[:,i]=rlat[:] 
    gradlat_2d[:,i]=gradlat[:] 
    
  lon_2d = np.ones(shape2d)
  gradlon_2d = np.ones(shape2d)
  for i in range(0,np.size(rlat)):
    lon_2d[i,:]=rlon[:]
    gradlon_2d[i,:]=gradlon[:]
    
  PHIt=lat_2d*radian  

# This has been edited from the original - typo, switch gradlat/lon
# as lat changes the metres change in proportion to the cosine
  e1t = radius * gradlon_2d * radian* np.cos(PHIt) 
  e2t = radius * gradlat_2d * radian 

# ==============================
# ==== u, v and f variants. ====
# u = i+0.5, j
# v = i, j+0.5
# f = i+0.5, j+0.5

  latv = rlat + 0.5*gradlat
  lonu = rlon + 0.5*gradlon 

  gradlatv = np.gradient(latv)
  gradlonu = np.gradient(lonu)

# ==== v points ====
# Same lon, different lat.
# However, e1 depends on lat! 
  for i in range(0,np.size(rlon)):
    lat_2d[:,i]=latv[:] 
    gradlat_2d[:,i]=gradlatv[:] 

  PHIv=lat_2d*radian  
    
  e1v = radius * gradlon_2d * radian* np.cos(PHIv) 
  e2v = radius * gradlat_2d * radian 

# ==== u points ====
# Same lat, different lon
  for i in range(0,np.size(rlat)):
    lon_2d[i,:]=lonu[:]
    gradlon_2d[i,:]=gradlonu[:]

  PHIt=lat_2d*radian  
    
  e1u = radius * gradlon_2d * radian* np.cos(PHIt) 
  e2u = e2t # i.e. e2u = radius * gradlat_2d * radian

# ==== f points ====
# different lat and lon (lonu and latv)

  e1f = radius * gradlon_2d * radian* np.cos(PHIv) 
  e2f = e2v # i.e. e2f = radius * gradlat_2d * radian 


# ===============================
# ==== unrotated coordinates ====
  rotated_cs = iris.coord_systems.RotatedGeogCS(37.5, 177.5)

  rlat_2d=np.ones(shape2d)
  for i in range(0,np.size(rlon)):
    rlat_2d[:,i]=rlat[:] 

  vlat_2d=np.ones(shape2d)
  for i in range(0,np.size(rlon)):
    vlat_2d[:,i]=latv[:] 
    
  rlon_2d=np.ones(shape2d)
  for i in range(0,np.size(rlat)):
    rlon_2d[i,:]=rlon[:]

  ulon_2d=np.ones(shape2d)
  for i in range(0,np.size(rlat)):
    ulon_2d[i,:]=lonu[:]

  glont, glatt = iris.analysis.cartography.unrotate_pole(
                                    rlon_2d, rlat_2d,
                                    rotated_cs.grid_north_pole_longitude,
                                    rotated_cs.grid_north_pole_latitude)

  glonu, glatu = iris.analysis.cartography.unrotate_pole(
                                    ulon_2d, rlat_2d,
                                    rotated_cs.grid_north_pole_longitude,
                                    rotated_cs.grid_north_pole_latitude)

  glonv, glatv = iris.analysis.cartography.unrotate_pole(
                                    rlon_2d, vlat_2d,
                                    rotated_cs.grid_north_pole_longitude,
                                    rotated_cs.grid_north_pole_latitude)

  glonf, glatf = iris.analysis.cartography.unrotate_pole(
                                    ulon_2d, vlat_2d,
                                    rotated_cs.grid_north_pole_longitude,
                                    rotated_cs.grid_north_pole_latitude)

# ==== WRITE TO FILE ====
  print ("Storing NEMO coordinates files to {}".format(outdir))
  ncfile = Dataset('{}/expand_amm15.coordinates.nc'.format(outdir),'w')

  y_dim = ncfile.createDimension('lon', np.size(rlon))     # Y
  x_dim = ncfile.createDimension('lat', np.size(rlat))     # X
  latr = ncfile.createVariable('lat', 'f4', ('lat',))
  lonr = ncfile.createVariable('lon', 'f4', ('lon',))
  latout = ncfile.createVariable('gphit', 'f4', ('lat','lon',))
  lonout = ncfile.createVariable('glamt', 'f4', ('lat','lon',))
  latoutu = ncfile.createVariable('gphiu', 'f4', ('lat','lon',))
  lonoutu = ncfile.createVariable('glamu', 'f4', ('lat','lon',))
  latoutv = ncfile.createVariable('gphiv', 'f4', ('lat','lon',))
  lonoutv = ncfile.createVariable('glamv', 'f4', ('lat','lon',))
  latoutf = ncfile.createVariable('gphif', 'f4', ('lat','lon',))
  lonoutf = ncfile.createVariable('glamf', 'f4', ('lat','lon',))

  e1out = ncfile.createVariable('e1t', 'f4', ('lat','lon',))
  e2out = ncfile.createVariable('e2t', 'f4', ('lat','lon',))
  e1outu = ncfile.createVariable('e1u', 'f4', ('lat','lon',))
  e2outu = ncfile.createVariable('e2u', 'f4', ('lat','lon',))
  e1outv = ncfile.createVariable('e1v', 'f4', ('lat','lon',))
  e2outv = ncfile.createVariable('e2v', 'f4', ('lat','lon',))
  e1outf = ncfile.createVariable('e1f', 'f4', ('lat','lon',))
  e2outf = ncfile.createVariable('e2f', 'f4', ('lat','lon',))

  latr[:] = rlat
  lonr[:] = rlon

  latout[:] = glatt
  lonout[:] = glont
  latoutu[:] = glatu
  lonoutu[:] = glonu
  latoutv[:] = glatv
  lonoutv[:] = glonv
  latoutf[:] = glatf
  lonoutf[:] = glonf

  e1out[:] = e1t
  e2out[:] = e2t
  e1outu[:] = e1u
  e2outu[:] = e2u
  e1outv[:] = e1v
  e2outv[:] = e2v
  e1outf[:] = e1f
  e2outf[:] = e2f

  ncfile.description = 'Expanded AMM15 coordinates file created for AMM15 to be pre smoothed by shapiro filter'

  now = datetime.now()
  current_time = now.strftime("%Y/%M/%d %H:%M:%S")
  repos = subprocess.run(['git', 'config', '--get', 'remote.origin.url'],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
  repos = repos.stdout.decode('utf-8').strip('\n')
  branch = subprocess.run(['git', 'rev-parse', '--abbrev-ref',  'HEAD'],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
  branch = branch.stdout.decode('utf-8').strip('\n')


  parser = argparse.ArgumentParser()
  script = parser.prog

  ncfile.history = "Created by {} from branch {} of {} on {} ".format(script,branch[:],repos[:],current_time)
  ncfile.pyversion = platform.python_version()
  ncfile.System = platform.system()
  ncfile.Release = platform.release()

  ncfile.close()




