
from __future__ import print_function
import scipy.ndimage
from astropy.io import fits 
import numpy
import megaradrp.datamodel as dm 
from .getspaxdim import getspaxdim

def convert(infile,arcsec_per_pixel=0.2,sigma_conv=8.,expansion_factor=5,writeout=None,overwrite=False,):
    """ 
    Convert a MEGARA fits RSS file format 
    to a traditional IFU fits data cube (3-d data array)
    
    
    infile: string with the rss fits filename
    
    
    arcsec_per_pixel: float, desired arcsec per pixel. It will find the closest arcsec_per_pixel of your input, but not exact
    sigma_conv: float. By default is 8. This parameter is the sigma of the gaussian convolution used to populate the square-pixels around 
    the center of each hexagon. Modify it by your own risk!
    expansion_factor: integer. By default is 5. The expansion of NAXIS2 to decrease the effect of irrational factor of the xy-sampling. 
    Modify it by your own risk. sigma_conv~sqrt(3)*expansion_factor
    writeout: string. filename of the fits output. By default the output fits is not written. 
    overwrite: bool. 
    
    usage -> 
    ----------------
    import megararss2cube
    cube=megararss2cube.convert('final_rss.fits')
    --------------------
    returns a fits HDUlist 
    
    --------------------------------------------
    TODO Include astrometric information such as angle and RA, DEC positions
    TODO is the flux well calibrated after all the transformations?
    ---------------------------------------------
    """
    
    PLATESCALE = 1.2120  # arcsec / mm
    rss = fits.open( infile )
    phdr = rss[1].header
    dhdr = rss[0].header
    data = rss[0].data
    
    conff=dm.read_fibers_extension(phdr)
    bundles_values=conff.bundles.keys()
    sky_bundles=[]
    for bundlei in bundles_values:
        if phdr["BUN%03d_T" % bundlei]=='SKY':
            sky_bundles.append(bundlei)
    
    w0 = dhdr['CRVAL1']  # reference wavelength
    try : dw = dhdr['CRDELT1']  # wavelength step
    except : dw = dhdr['CDELT1']  # wavelength step
    wunit = dhdr['CUNIT1']  # wavelength unit
    wtype = 'WAVE'  # type spectra

    # define the dimensions of the spaxel array 
    Nx, Ny, x0, y0, dx, dy  = getspaxdim( data,phdr,sky_bundles,expansion_factor=expansion_factor)

    nbin=int(round(float(arcsec_per_pixel)/float(dx)))


    Nw = dhdr['NAXIS1']    # number of wave. steps
    

       
    # initialize an empty 3-d cube (zero everywhere)
    cube = fits.PrimaryHDU()
    cube.header.update(NAXIS=3)
    cube.header.update(NAXIS1=Nx)
    cube.header.update(NAXIS2=Ny)
    cube.header.update(NAXIS3=Nw)
    cube.header.update(CD1_1=dx)
    cube.header.update(CD2_2=dy)
    cube.header.update(CD3_3=dw)
    cube.header.update(CRPIX1=0)
    cube.header.update(CRPIX2=0)
    cube.header.update(CRPIX3=0)
    cube.header.update(CRVAL1=x0)
    cube.header.update(CRVAL2=y0)
    cube.header.update(CRVAL3=w0)

    cube.header.update(CTYPE1='RA---CAR')
    cube.header.update(CTYPE2='DEC--CAR')
    cube.header.update(CTYPE3=wtype)
    cube.header.update(CUNIT3=wunit)

    cube.header.update(CD1_2=0)
    cube.header.update(CD1_3=0)
    cube.header.update(CD2_1=0)
    cube.header.update(CD2_3=0)
    cube.header.update(CD3_1=0)
    cube.header.update(CD3_2=0)


    cube.data = numpy.zeros( (Nw,Ny,Nx) )

    # extract each spectrum and place it
    # into the 3-d cube
    for ispec in range(len(data)): 
     fib_str='{:3d}'.format(ispec+1)
     fib_str=fib_str.replace(' ','0') 
     if not(phdr['FIB'+fib_str+'_B'] in sky_bundles):
        try:
         end_sp=phdr['FIB'+fib_str+'W2'] 
         start_sp=phdr['FIB'+fib_str+'W1']
        except:
          if ('start_sp' in locals()):
           print('Warning! FIB'+fib_str+'W1 and W2 information missing in header. Assuming previous fiber wavelength coverage.') 
          else: 
           end_sp=Nw
           start_sp=1   
           print('Warning! FIB'+fib_str+'W1 and W2 information missing in header. Assuming default wavelength coverage.') 

        spec = data[ispec][:]
        Nwspec = Nw    

        xpos = (phdr['FIB'+fib_str+'_x']+5.)*PLATESCALE  
        ypos = (phdr['FIB'+fib_str+'_y']+5.)*PLATESCALE
        ix = int( round((xpos - x0),3) / dx  )
        iy = int( round((ypos - y0),3) / dy  )

        lambda_arr=w0+dw*numpy.arange(0,Nwspec,1)

        for i in range( start_sp, min(end_sp,Nwspec) ):

            

            cube.data[i][iy][ix] = spec[i]*3.00e-5/lambda_arr[i]**2 ## Jy to erg/s/cm**2/A
    print('1st step')   
    for i in range( start_sp, min(end_sp,Nwspec) ):
        print(str(i)+'/'+str(Nwspec)+' spectral channels',end="\r")
        cube.data[i]=scipy.ndimage.filters.gaussian_filter(cube.data[i], sigma=sigma_conv)
    
    
    cube_rebin = fits.PrimaryHDU()
    cube_rebin.header.update(NAXIS=3)
    cube_rebin.header.update(NAXIS1=Nx//nbin)
    cube_rebin.header.update(NAXIS2=Ny//nbin)
    cube_rebin.header.update(NAXIS3=Nw)
    cube_rebin.header.update(CD1_1=dx*nbin)
    cube_rebin.header.update(CD2_2=dy*nbin)
    cube_rebin.header.update(CD3_3=dw)
    cube_rebin.header.update(CRPIX1=0)
    cube_rebin.header.update(CRPIX2=0)
    cube_rebin.header.update(CRPIX3=0)
    cube_rebin.header.update(CRVAL1=x0)
    cube_rebin.header.update(CRVAL2=y0)
    cube_rebin.header.update(CRVAL3=w0)
        
    cube_rebin.header.update(CTYPE1='RA---CAR')
    cube_rebin.header.update(CTYPE2='DEC--CAR')
    cube_rebin.header.update(CTYPE3=wtype)
    cube_rebin.header.update(CUNIT3=wunit)
        
    cube_rebin.header.update(CD1_2=0)
    cube_rebin.header.update(CD1_3=0)
    cube_rebin.header.update(CD2_1=0)
    cube_rebin.header.update(CD2_3=0)
    cube_rebin.header.update(CD3_1=0)
    cube_rebin.header.update(CD3_2=0)
    cube_rebin.header.update(BUNIT= 'erg/s/cm**2/Angstrom')                                      



        
    cube_rebin.data = numpy.zeros( (Nw,Ny//nbin,Nx//nbin) )
    print('')
    print('2nd step')
    for i in range( 0, Nwspec) :        
     shape=cube.data[i].shape    
     print(str(i)+'/'+str(Nwspec)+' spectral channels',end="\r")
     for xi in numpy.arange(0,shape[0],nbin)[:-1]:
         for yj in numpy.arange(0,shape[1],nbin)[:-1]:
             pixel_ij=numpy.sum(cube.data[i][xi:xi+nbin,yj:yj+nbin])             
             cube_rebin.data[i][xi//nbin,yj//nbin]=pixel_ij     
    if writeout !=None:
        cube_rebin.writeto(writeout,overwrite=overwrite)
    return( cube_rebin)
