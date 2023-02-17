# Surge data description

Surge LAT data comes in txt format from Colin Bell 

Two txt for for CS2X and CS20

## Format of data

Script to injest ascii data from ColinBell of LAT 
Note format is :

```
2  2 62.83333 -19.75000  1.51 -1.51  2.69  2.02  1.00  0.33



values 1,2: internal coordinate
values 3,4: latitude, longitude
value 5: HAT
value 6: LAT (usually negative is it is below MSL)
values 7-10: MHWS, MHWN, MLWN, MLWS�|  (all the LAT datum)
```


So we can create variables for all values given test if the lat,lons are truly rectangular



