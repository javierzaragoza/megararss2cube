# megararss2cube

Megararss2cube is a tool to convert MEGARA LCB reduced dataproducts from the RSS format obtained with [megaradrp](https://github.com/guaix-ucm/megaradrp) to a more user-friendly 3d datacube. 

### Install

Download the [source code](https://github.com/javierzaragoza/megararss2cube/archive/master.zip) and extract it. Then:

    python setup.py install


### Dependencies

scipy
astropy
numpy
megaradrp

### Usage

    import megararss2cube
    #convert a rss fits dataproduct 'rss_product.fits' to a cube hdulist
    #optionally, if writeout is defined, it writes the cube product to 'cube_product.fits'
    cube=megararss2cube.convert('rss_product.fits',writeout='cube_product.fits')

### Contact

javier.zaragoza@inaoep.mx
