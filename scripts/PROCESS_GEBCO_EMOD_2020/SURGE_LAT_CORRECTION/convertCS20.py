"""
Script to injest ascii data from ColinBell of LAT 
Note format is :
2  2 62.83333 -19.75000  1.51 -1.51  2.69  2.02  1.00  0.33

values 1,2: internal coordinate
values 3,4: latitude, longitude
value 5: HAT
value 6: LAT (usually negative is it is below MSL)
values 7-10: MHWS, MHWN, MLWN, MLWS  (all the LAT datum)


So we can create variables for all values given test if the lat,lons are truly rectangular
"""

from netCDF4 import Dataset   # for writing out the data
import numpy    as np
import numpy.ma as ma

numlines=sum(1 for line in open('cs20_stat.txt')) # a
print (numlines)


fin=open('cs20_stat.txt')
data = (np.fromfile(fin,sep=" ",count=numlines*10,dtype=np.float))
data= data.reshape(numlines,10)
y         = (data[:,1])
x         = (data[:,0])
print(np.max(y),np.max(x),np.min(x),np.min(y))

X=1+int(np.max(x)-np.min(x))  # range doesn't start at zero so account for that using both max and min
Y=1+int(np.max(y)-np.min(y))

y=y-np.min(y) # start at 0
x=x-np.min(x) # start at 0

latitude  = np.ones( (Y,X) ) * -999.
longitude = np.ones( (Y,X) ) * -999.
HAT       = np.ones( (Y,X) ) * -999.
LAT       = np.ones( (Y,X) ) * -999.
#Want alternative mask (+999) for land poiints intrewad of poitns in the saea but missing
LAT[:,456:]      = +999
LAT[:642,192:]      = +999
LAT[512:630,162:211]      = +999
LAT[189:450,110:216]      = +999
LAT[205:355,45:120]      = +999

LAT[351:381,61:116]      = +999

LAT[708:772,406:452]      = +999
LAT[727:732,388:399]      = +999
LAT[627:687,321:391]      = +999
MHWS      = np.ones( (Y,X) ) * -999.
MHWN      = np.ones( (Y,X) ) * -999.
MLWN      = np.ones( (Y,X) ) * -999.
MLWS      = np.ones( (Y,X) ) * -999.


for line in range(0,numlines):

#Looks like indexes start at top, so reverse here for ease of viewing in ncview etc.

  latitude  [ (Y-1)-int(y[line]), int(x[line]) ]  = data[line,2]
  longitude [ (Y-1)-int(y[line]), int(x[line]) ]  = data[line,3]
  HAT       [ (Y-1)-int(y[line]), int(x[line]) ]  = data[line,4]
  LAT       [ (Y-1)-int(y[line]), int(x[line]) ]  = data[line,5]
  MHWS      [ (Y-1)-int(y[line]), int(x[line]) ]  = data[line,6]
  MHWN      [ (Y-1)-int(y[line]), int(x[line]) ]  = data[line,7]
  MLWN      [ (Y-1)-int(y[line]), int(x[line]) ]  = data[line,8]
  MLWS      [ (Y-1)-int(y[line]), int(x[line]) ]  = data[line,9]


## Looks like a regular square grid
## se we need not mask out lat lons
#just take the max over row/column
for i in range(0,Y):
 for j in range(0,X):
    latitude  [i,j] = np.max (latitude [i,:])
    longitude [i,j] = np.max (longitude[:,j])


### Now write to file

root_grp = Dataset('cs20_stat.nc', 'w')

## dimensions

root_grp.createDimension('x', X )
root_grp.createDimension('y', Y )

## variables

latitude_OUT  = root_grp.createVariable('nav_lat','f4',('y','x',))
longitude_OUT = root_grp.createVariable('nav_lon','f4',('y','x',))
HAT_OUT       = root_grp.createVariable('HAT'    ,'f4',('y','x',))
LAT_OUT       = root_grp.createVariable('LAT'    ,'f4',('y','x',))
MHWS_OUT      = root_grp.createVariable('MHWS'   ,'f4',('y','x',))
MHWN_OUT      = root_grp.createVariable('MHWN'   ,'f4',('y','x',))
MLWN_OUT      = root_grp.createVariable('MLWN'   ,'f4',('y','x',))
MLWS_OUT      = root_grp.createVariable('MLWS'   ,'f4',('y','x',))

##
LAT[0,:] = -999
LAT[:,0] = -999
LAT[-1,:] = -999
#LAT[:,0:10] = -999

latitude_OUT [:]  = latitude [:]
longitude_OUT[:]  = longitude[:]
HAT_OUT      [:]  = ma.masked_equal(HAT [:],-999)

LAT_OUT      [:]  = ma.masked_equal(LAT [:],-999)
MHWS_OUT     [:]  = ma.masked_equal(MHWS[:],-999)
MHWN_OUT     [:]  = ma.masked_equal(MHWN[:],-999)
MLWN_OUT     [:]  = ma.masked_equal(MLWN[:],-999)
MLWS_OUT     [:]  = ma.masked_equal(MLWS[:],-999)

root_grp.close()
## Want outer permineter masked as want no extrapolation beyond it later
