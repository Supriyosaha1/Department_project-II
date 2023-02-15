## Simple routines to write parameter files
import numpy as np

# ----------------------------------------------------------------------
# functions to write RASCAS parameter files from dictionaries.
# ----------------------------------------------------------------------

class param_section():
    def __init__(self,name,params):
        self.name = name
        self.params = params


def write_param_section(f,s):
    f.write('[%s] \n'%s.name)
    for k in s.params.keys():
        f.write('  %s = %s \n'%(k,s.params[k]))
    f.write('\n') # skip a line between sections
        
def write_parameter_file(filename,params):
    f = open(filename,'w')
    for s in params:
        write_param_section(f,s)
    f.close()
    
def print_param_section(ps):
    print('[%s] \n'%ps.name)
    for k in ps.params.keys():
        print('  %s = %s \n'%(k,ps.params[k]))
        
        
# -------------------------------------------------------------------------
# functions to generate parameters for mock observations using peeling off. 
# -------------------------------------------------------------------------
class fluxParams():
    def __init__(self,aperture):
        self.aperture = aperture  # this is radius of circular aperture within which we collect photon packets. [code units]
    def write(self,f):
        f.write("%.10e \n "%(self.aperture))
    def toString(self):
        return "%.10e \n "%(self.aperture)
    
class specParams():
    def __init__(self,npix,aperture,lmin,lmax):
        self.npix     = npix
        self.aperture = aperture
        self.lmin     = lmin
        self.lmax     = lmax
    def write(self,f):
        f.write("%i %.10e %.10e %.10e \n "%(self.npix,self.aperture,self.lmin,self.lmax))
    def toString(self):
        return "%i %.10e %.10e %.10e \n "%(self.npix,self.aperture,self.lmin,self.lmax)
        
class imageParams():
    def __init__(self,npix,side):
        self.npix = npix
        self.side = side
    def write(self,f):
        f.write("%i %.10e \n "%(self.npix,self.side))
    def toString(self):
        return "%i %.10e \n "%(self.npix,self.side)

class cubeParams():
    def __init__(self,cube_lbda_npix,cube_image_npix,cube_lmin,cube_lmax,cube_side):
        self.cube_lbda_npix  = cube_lbda_npix  
        self.cube_image_npix = cube_image_npix 
        self.cube_lmin       = cube_lmin       
        self.cube_lmax       = cube_lmax       
        self.cube_side       = cube_side       
    def write(self,f):
        f.write("%i %i %.10e %.10e %.10e \n "%(self.cube_lbda_npix,self.cube_image_npix,self.cube_lmin,self.cube_lmax,self.cube_side))
    def toString(self):
        return "%i %i %.10e %.10e %.10e \n "%(self.cube_lbda_npix,self.cube_image_npix,self.cube_lmin,self.cube_lmax,self.cube_side)

class pointing():
    # a pointing is a direction of observations, along which we'll possibly collect a flux, a spectrum, an image, and/or a cube. 
    def __init__(self,k,c,flux=None,spec=None,image=None,cube=None):
        self.kobs = np.array(k) / np.sqrt(k[0]**2+k[1]**2+k[2]**2) # make sure direction vector is normalised to 1
        self.targetpos = c
        # define photometric observation 
        if flux is not None:
            self.flux = flux
        else:
            self.flux = fluxParams(0)
        # define spectroscopic observation
        if spec is not None:
            self.spec = spec
        else:
            self.spec = specParams(0,0,0,0)
        # define imaging observation
        if image is not None:
            self.image = image
        else:
            self.image = imageParams(0,0)
        # define 3D spectroscopy observation
        if cube is not None:
            self.cube = cube
        else:
            self.cube = cubeParams(0,0,0,0,0)
        
    def write(self,f):
        f.write("%.10e %.10e %.10e \n "%(self.kobs[0],self.kobs[1],self.kobs[2]))
        f.write("%.10e %.10e %.10e \n "%(self.targetpos[0],self.targetpos[1],self.targetpos[2]))
        self.flux.write(f)
        self.spec.write(f)
        self.image.write(f)
        self.cube.write(f)
        
    def toString(self):
        s = "%.10e %.10e %.10e \n "%(self.kobs[0],self.kobs[1],self.kobs[2])
        s = s + "%.10e %.10e %.10e \n "%(self.targetpos[0],self.targetpos[1],self.targetpos[2])
        s = s + self.flux.toString()
        s = s + self.spec.toString()
        s = s + self.image.toString()
        s = s + self.cube.toString()
        return s
        
