# Script to include useful tools for processing LSM. 
# Including:
# remove_lakes - find and fill isolated ocean points.
# find_1pt_islands - find 1 pt islands.
# find_1pt_estuaries - find 1 pt estuaries.
# flip - replace land/sea convension, change 1 = 0 and 0 = 1.
# fill - replace missing data with nearest neighbour.
#
# NB. if any changes made, make sure delete .pyc file.
#
# JGraham - created 2-Sep-2015

################################
# Remove lakes from the mask
################################

def remove_lakes(maskin):
    """
    Remove isolated 'lakes' from a land/sea mask

    INPUT = maskin, NB. land = 0, sea = 1.

    OUTPUT = maskout, with lakes removed. 
    """

    import numpy as np

    # Sweep grid and identify each individual basin.
    max_nbr_of_basins = 500 
    # memory error in python if >580 on local machine
    print ('Max number of basins: ', max_nbr_of_basins)

    maskin[0,:] = 0
    maskin[:,0] = 0
    maskin[-1,:] = 0
    maskin[:,-1] = 0

    [lm,mm] = np.shape(maskin)
    lm=lm-2
    mm=mm-2

    tmp = np.zeros((lm + 4, mm + 4))
    tmp[1:-1, 1:-1] = maskin
    maskin = tmp 

    ## do not forget to remove margins at the end! ##

    wet = np.zeros(np.shape(maskin))
    Relations = np.zeros((max_nbr_of_basins,lm*mm))
    Relations1= np.zeros((max_nbr_of_basins,lm*mm)) # needed if >max basins
    cells_to_add = np.zeros((lm, mm))

    # Initialization.
    NoBasin  = 0
    NoBasin1 = 0

    wet[maskin > 0] = 1
    n=0  
    for j in range(mm): # inner lon
       if j%100==0:
          print (j)
       for i in range(lm): # inner lat
          # Must use ic,jc in place of i,j (inner domain).
          ic = i + 2
          jc = j + 2

          # If a cell is wet and if it is not related to a basin already identified,
          #  then you just discovered a new basin.
          if ( wet[ic, jc] and (~np.any(Relations[:, lm*j+i])) 
                  and  (~np.any(Relations1[:, lm*j+i]))    
                  and  (NoBasin < max_nbr_of_basins)):
             NoBasin = NoBasin + 1 

             if (NoBasin == max_nbr_of_basins):
                print ('Too many basins, move to overflow array...')
                Relations1 = np.copy(Relations)
                Relations = np.zeros((max_nbr_of_basins,lm*mm))
                NoBasin1 = NoBasin
                NoBasin = 1
                n = n+1
                print ('n = ', n)

             if (NoBasin < max_nbr_of_basins):
                if NoBasin%100==0:
                   print (NoBasin, ' basins... ')

                # The cell (i,j) must be added to the basin # NoBasin.
                cells_to_add[:]    = 0
                cells_to_add[i, j] = 1 

                # So you have discovered a new basin. Now, you must define the extent
                #  by sweeping the grid and adding all cells that are part of this basin.

                not_done = 1
                while not_done: 
                   
                   # If no more cells to add/basin fully explored, look for a new basin.
                   if ( ~np.any(cells_to_add) ):
                      not_done = 0
                      break

                   # inner_j_loop: for ix2 = inner lon dimension
                   for ix2 in range(mm):
 
                      if np.any( cells_to_add[:, ix2] ):  # To accelerate sweep
 
                         # inner_i_loop: for ix1 = 1, lm
                         for ix1 in range(lm):

                            if ( cells_to_add[ix1, ix2] ):

                               Relations[NoBasin, lm*ix2+ix1] = 1
                               cells_to_add[ix1, ix2] = 0

                               # Look for adjacent wet cells that are not related elsewhere
                               # NB. wet array covers whole domain, not just inner domain (+2)
                               if (wet[ix1 + 3, ix2 + 2]):
                                  if (Relations[NoBasin, lm*ix2 + (ix1+1)]==0):
                                     cells_to_add[ix1 + 1, ix2] = 1

                               if (wet[ix1 + 1, ix2 + 2]):
                                  if (Relations[NoBasin, lm *ix2 + (ix1 - 1)]==0):
                                     cells_to_add[ix1 - 1, ix2 ] = 1

                               if ( wet[ix1 + 2, ix2 + 3] ):
                                  if (Relations[NoBasin, lm *(ix2 + 1) + ix1]==0):
                                     cells_to_add[ix1    , ix2 + 1] = 1

                               if ( wet[ix1 + 2, ix2 + 1] ):
                                  if (Relations[NoBasin, lm*(ix2 - 1) + ix1]==0):
                                     cells_to_add[ix1    , ix2 - 1] = 1


    # Find the main basin (the one with the largest area).
    print ('Finished exploring the grid...'  )
    AreaMainBasin = 0
 
    for i in range(max_nbr_of_basins):
       if (np.sum(Relations[i,:]) > AreaMainBasin ):
          AreaMainBasin = np.sum(Relations[i,:])
       if NoBasin1:
          if (np.sum(Relations1[i,:]) > AreaMainBasin ):
             AreaMainBasin = np.sum(Relations1[i,:])      

    if NoBasin1:
       NoBasin=NoBasin+NoBasin1   

    print ('There are ', NoBasin, ' different basins over the grid.')
    print ('Area of main basin = ', AreaMainBasin, ' cells.'       )
  
    # Consider all other basins as `lakes', and fill them.
    for i in range(max_nbr_of_basins):
       if (np.sum(Relations[i, :]) < AreaMainBasin):
          j = np.squeeze(np.where( Relations[i, :] == 1 ))
          if ( np.size(j) > 0 ):
             i_h = j - lm * np.floor(j/ lm)
             j_h = np.floor(j/lm)
             maskin[(i_h.astype(int) + 2), (j_h.astype(int) + 2)] = 0

    # Check additional array if max basins exceeded   
    if NoBasin1:
       for i in range(max_nbr_of_basins):
          if (np.sum(Relations1[i, :]) < AreaMainBasin):
             j = np.squeeze(np.where( Relations1[i, :] == 1 ))
             if ( np.size(j) > 0 ):
                i_h = j - lm * np.floor(j/ lm)
                j_h = np.floor(j/lm)
                maskin[(i_h.astype(int) + 2), (j_h.astype(int) + 2)] = 0
      
    # Remove bry margins
    maskout = maskin[1:-1,1:-1]
    maskout[0,:] = maskout[1,:]
    maskout[:,0] = maskout[:,1]
    maskout[-1,:] = maskout[-2,:]
    maskout[:,-1] = maskout[:,-2]
    return(maskout)   


################################
# Find 1 pt. islands
################################
def find_1pt_islands(maskin):
    """
    Find 1 pt islands, completely separate from land.
    Ignore any diagonal connections.
    
    INPUT = maskin, NB. land = 1, sea = 0
    OUTPUT = flag (island = 1)
    """
    import numpy as np
    from pylab import nan
 
    nlat,nlon = np.shape(maskin) 
    area=np.ones((nlat,nlon))*nan
    flag=np.zeros((nlat,nlon))

    n=0
    for i in range(1,nlat-1):
       if i%100==0:
          print (i, ' of ', nlat)

       for j in range(1,nlon-1):

          # if land, calculate area of surrounding land
          if maskin[i,j]:
             area[i,j] = (maskin[i-1,j]+maskin[i+1,j]+
                           maskin[i,j+1]+maskin[i,j-1]+
                           maskin[i-1,j-1]+maskin[i-1,j+1]+
                           maskin[i+1,j-1]+maskin[i+1,j+1])

             if area[i,j]==0:
                n=n+1
                flag[i,j]=1

    print ('total number of 1 pt. islands: ', n)
    return(flag)

################################
# Find 1 pt. estuaries
################################
def find_1pt_estuaries(maskin):
    """
    Find 1 pt estuaries, including diagonal connections.
    
    INPUT = maskin, NB. land = 1, sea = 0
    OUTPUT = flag (estuary pt. = 1)

    Inlcudes a quick fix to also flag 1pt lakes produced.
    NB. will likely create more lakes -> remove_lakes again.
    """
    import numpy as np
    from pylab import nan
 
    nlat,nlon=np.shape(maskin)

    widthi=np.ones((nlat,nlon))*nan
    widthj=np.ones((nlat,nlon))*nan

    widthi[1:-1,:]=maskin[0:-2,:]+maskin[2:,:]
    widthj[:,1:-1]=maskin[:,0:-2]+maskin[:,2:]

    # mask points = land
    widthi[maskin==1] = 0
    widthj[maskin==1] = 0

    flag=np.zeros((nlat,nlon))

    #####
    # Try a simple method, just patching all narrow channels (land on either side)
    #####
    flag[widthi==2] = 1
    flag[widthj==2] = 1

    ######
    # This section looks for diagonal estuaries - enclosed on two sides, 
    # but not just in one direction across each pt.
    # e.g. for 1 = land, 0 = sea 
    #
    # 1 0 0 0 0 1 1
    # 1 1 0 0 1 1 1
    # 1 1 1 0 0 1 1
    # 1 1 1 1 0 0 1
    #
    # Loop over domain and check for water access to each point?
    #
    ######
    access = 4-(widthi+widthj)
    access[maskin==1] = 0

    for i in range(1,nlat-1):
       if i%100==0:
          print (i, ' of ', nlat)
       for j in range(1,nlon-1):

          # Check access to point (if not already flagged)
          if access[i,j]==2 and flag[i,j]==0:

             # check for adjacent ocean points
             if (maskin[i-1,j]==0 and maskin[i,j-1]==0 
                and maskin[i-1,j-1]==1): 
                flag[i,j]=1

             elif (maskin[i-1,j]==0 and maskin[i,j+1]==0 
                and maskin[i-1,j+1]==1): 
                flag[i,j]=1
            
             elif (maskin[i+1,j]==0 and maskin[i,j+1]==0 
                and maskin[i+1,j+1]==1): 
                flag[i,j]=1

             elif (maskin[i+1,j]==0 and maskin[i,j-1]==0 
                and maskin[i+1,j-1]==1): 
                flag[i,j]=1

    # Quick fix to also flag one point lakes that have been created?
    ### Note, will likely be left with other lakes -> run remove_lakes again. ###
    mask2 = np.copy(maskin)
    mask2[flag==1] = 1

    widthi = np.ones((nlat,nlon))*nan
    widthj = np.ones((nlat,nlon))*nan

    widthi[1:-1,:] = mask2[0:-2,:]+mask2[2:,:]
    widthj[:,1:-1] = mask2[:,0:-2]+mask2[:,2:]

    area = widthi+widthj
    area[mask2==1] = 0

    flag[area==4] = 1

    # Remove points in contact with the wider ocean 
    # (i.e. keep 1 point inlets, and
    # longer estuaries will be left as 1 pt inlets)
    mask2 = np.copy(maskin)
    mask2[flag==1] = 1    
    area = np.ones((nlat,nlon))*nan 
    area[1:-1,1:-1] = (mask2[0:-2,1:-1]+mask2[2:,1:-1]+
                     mask2[1:-1,0:-2]+mask2[1:-1,2:])        

    flag[area<4] = 0

    total = np.sum(flag[:])
    print ('total number of estuary points: ', total)
    
    return(flag)


################################
# Flip mask, so 1 = 0 and 0 = 1
################################
def flip(maskin):
    """
    Flip the mask format, e.g.:
    
    INPUT = maskin, with land = 1, sea = 0
    OUTPUT = maskout, with land = 0, sea = 1
    (or vice versa)
    """
    import numpy as np
    from pylab import nan

    maskout = np.copy(maskin)
    maskout[maskin==1]=nan
    maskout[maskin==0]=1
    maskout[np.isnan(maskout)]=0

    return(maskout)


################################
# Fill masked data 
################################
def fill(data, invalid=None):
    """
    From Enda: in EXTEND_SPLICE_AND_FILL_TPXO.py
    Replace the value of invalid 'data' cells (indicated by 'invalid')
    by the value of the nearest valid data cell
    Input:
        data:    numpy array of any dimension
        invalid: a binary array of same shape as 'data'. True cells set where data
                 value should be replaced.
                 If None (default), use: invalid  = np.isnan(maskin)
    Output:
        Return a filled array.
    """    

    import numpy as np
    import scipy.ndimage as nd
    if invalid is None: 
        invalid = np.isnan(data)

    ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    filled = data[tuple(ind)]
    return filled


################################
# rx0 factor? 
################################
def rfactor(h, mask=None):
    """
    Adapted from ROMS function. 
    Input:
        h:       bathymetry
        invalid: Land-Sea Mask (expect 1 = land)
                 If None (default), use: mask  = np.isnan(h)
    Output:
        Return an array of rx0 values, same size as h.
    """    

    import numpy as np
    from pylab import nan

    if mask is None: 
        mask = ~np.isnan(h)

    L,M = np.shape(h)

    umask = np.zeros([L,M])
    vmask = np.zeros([L,M])
    
    # Land/Sea mask on U-points.
    for j in range(M):
        for i in range(1,L):
            umask[i-1,j]=mask[i,j]*mask[i-1,j]

    # Land/Sea mask on V-points.
    for j in range(1,M):
        for i in range(L):
            vmask[i,j-1]=mask[i,j]*mask[i,j-1]

    #------------------------------------------------------------
    #  Compute R-factor.
    #------------------------------------------------------------
    hx = np.zeros([L,M])
    hy = np.zeros([L,M])

    tmp = h[1:,:]+h[:-1,:]
    tmp[tmp==0] = nan
    hx[:-1,:]=np.abs(h[1:,:]-h[:-1,:])/(tmp)

    tmp = h[:,1:]+h[:,:-1]
    tmp[tmp==0] = nan
    hy[:,:-1]=np.abs(h[:,1:]-h[:,:-1])/(tmp)
    
    hx[np.isnan(hx)] = 0
    hy[np.isnan(hy)] = 0
    hx=hx*umask
    hy=hy*vmask

    r=np.maximum(np.maximum(hx[:-1,:-1],hx[:-1,1:]), 
               np.maximum(hy[:-1,:-1],hy[1:,:-1]))

    rmin=np.min(r)
    rmax=np.max(r)
    ravg=np.nanmean(r)
    rmed=np.median(r)

    print ('Minimum r-value = ', str(rmin))
    print ('Maximum r-value = ', str(rmax))
    print ('Mean    r-value = ', str(ravg))
    print ('Median  r-value = ', str(rmed))

    return r


