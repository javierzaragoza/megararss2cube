import numpy

def getspaxdim( data ,phdr,sky_bundles,expansion_factor=5):
    """
    define the dimensions of the spaxel array. Return values are
      Nx : number of spaxels in the x direction
      Ny : number of spaxels in the y direction
      dx : arcseconds between adjacent spaxels in the x direction
      dy : arcseconds between adjacent spaxels in the y direction
    
    MEGARA has hexagonal spaxels. To convert them from hexagons 
    to squares a factor of numpy.sqrt(3) must be applied to NAXIS1. 
    Then, the spatial grid x,y needs to be expanded so the irrational
    factor does not affect the spatial sampling. By default is 5, 
    play with it by your own risk. The dimensions are then too large, 
    but this will be fixed rebinning the very large cube in the 
    2nd step.
    


    """
    PLATESCALE = 1.2120  # arcsec / mm MEGARA PLATESCALE
    # reduce a list to just the unique elements
    
    uniquelist = lambda x : dict(map(lambda i: (i,1),x)).keys()

    Nspec = len(data)

    xpos, ypos = [], []
    # collect all the physical spaxel locations
    # rounded to the nearest 0.001 arcseconds
    for ispec in range(Nspec): 
        fib_str='{:3d}'.format(ispec+1)
        fib_str=fib_str.replace(' ','0') 
 
        if not(phdr['FIB'+fib_str+'_B'] in sky_bundles):
         xpos.append(round((phdr['FIB'+fib_str+'_x']+5.)*PLATESCALE,3))    
         ypos.append(round((phdr['FIB'+fib_str+'_y']+5.)*PLATESCALE,3))
        


    # reduce to lists of just the unique x,y positions
    
    xpos = sorted( uniquelist( xpos ) )
    ypos = sorted( uniquelist( ypos ) )

    # find the separation between positions, 
    # rounded to the nearest 0.001 arcseconds
    xsteps = [ round(xpos[i]-xpos[i-1],3) for i in range(1,len(xpos)) ]
    ysteps = [ round(ypos[i]-ypos[i-1],3) for i in range(1,len(ypos)) ]

    # reduce to lists of unique separations
    xsteps = sorted( uniquelist( xsteps ) )
    ysteps = sorted( uniquelist( ysteps ) )

    
    # we probably have a single value in x and y now,
    # and for symmetric spaxels these are probably the same.
    # We take the minimum as defining the unit of separation
    # between adjacent spaxels

    dxspax = min( xsteps )
    dyspax = min( ysteps )
    
    # Finally, we assume the spaxel array spans the 
    # total x and y spatial extent with NXSPAX and NYSPAX
    # equally spaced spaxels.  Those dimensions are returned.
    NXSPAX = int( (max(xpos) - min(xpos))/dxspax ) + 1 
    NYSPAX = int( (max(ypos) - min(ypos))/dyspax ) + 1

    
    factor_Nxpix=numpy.sqrt(3)*expansion_factor
    return( int(round(NXSPAX*factor_Nxpix)), NYSPAX*expansion_factor, min(xpos), min(ypos), dxspax/factor_Nxpix, dyspax/expansion_factor)
