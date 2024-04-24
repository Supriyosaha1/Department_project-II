
# -*- coding: utf-8 -*-
import numpy as np
from scipy.io import FortranFile as ff 
import domain as d

class gas(object):
    """ This class defines gas objects.
    """
    
    def __init__(self,mix,nhi,dopwidth,vleaf,ndust,xleaf,leaflevel,vturb):
       
        self.mix = mix
        if self.mix == 'HI':
            self.nhi       = nhi
            self.dopwidth  = dopwidth
            self.vleaf     = vleaf
            self.xleaf     = xleaf
            self.leaflevel = leaflevel
        elif self.mix == ('HI_D_dust' or 'HI_dust'):
            self.nhi       = nhi
            self.dopwidth  = dopwidth
            self.vleaf     = vleaf
            self.ndust     = ndust
            self.xleaf     = xleaf
            self.leaflevel = leaflevel
        elif self.mix == 'SiII_dust':
            self.nSiII     = nhi
            self.dopwidth  = dopwidth
            self.vleaf     = vleaf
            self.ndust     = ndust
            self.xleaf     = xleaf
            self.leaflevel = leaflevel
        elif self.mix == 'new_ions':
            self.nion      = nhi
            self.dopwidth  = dopwidth
            self.vturb     = vturb
            self.vleaf     = vleaf
            self.ndust     = ndust
            self.xleaf     = xleaf
            self.leaflevel = leaflevel
        elif self.mix == 'dust':
            self.dopwidth  = dopwidth
            self.vturb     = vturb
            self.vleaf     = vleaf
            self.ndust     = ndust
            self.xleaf     = xleaf
            self.leaflevel = leaflevel
           
class mesh(object):
    """ This class manages mesh objects, which contain a domain object, the mesh data structure itself, and a gas object.
    """
    
    def __init__(self, filename=None, gasmix=None, Silent=True):

        if filename is None:
            self.domain   = None
            self.ncoarse  = 0
            self.noct     = 0
            self.nleaf    = 0
            self.father   = None
            self.son      = None
            self.nbor     = None
            self.octlevel = None
            self.xoct     = None
            self.xleaf    = None
            self.gas      = None
        else:
            if not(Silent): print("-----> reading mesh in file ",filename)
            f = ff(filename)
            
            # read domain
            #domainType = str(f.read_record('a10')[0])
            #domainType = domainType.strip()
            domainType = (f.read_record('S10')[0]).decode('UTF-8').strip()
            if not(Silent):
                print("-----> domain")
                print("domain type =", domainType)
            if domainType == 'sphere':
                center = f.read_reals('d')
                radius = f.read_reals('d')
                if not(Silent): print("domain size =",radius)
                self.domain    = d.sphere(center,radius)
            elif domainType == 'shell':
                center = f.read_reals('d')
                rin    = f.read_reals('d')
                rout   = f.read_reals('d')
                if not(Silent): print("domain size =",rin,rout)
                self.domain    = d.shell(center=center,rin=rin,rout=rout)
            elif domainType == 'cube':
                center = f.read_reals('d')
                size   = f.read_reals('d')
                self.domain    = d.cube(center,size)
            elif domainType == 'slab':
                zc = f.read_reals('d')
                thickness = f.read_reals('d')
                center = zc
                self.domain = d.slab(zc,thickness)
            else:
                return IOError("type not defined",domainType)
            if not(Silent): print("domain center =",center)
    
            # read mesh
            [self.ncoarse,self.noct,self.ncell,self.nleaf] = f.read_ints()
            if not(Silent): 
                print("-----> mesh")
                print("ncoarse =",self.ncoarse)
                print("noct    =",self.noct)
                print("ncell   =",self.ncell)
                print("nleaf   =",self.nleaf)
            # father
            self.father = f.read_ints()
            if not(Silent): print("INFO father   ",np.shape(self.father),np.min(self.father),np.max(self.father))
            # son
            self.son = f.read_ints()
            if not(Silent): print("INFO son      ",np.shape(self.son),np.min(self.son),np.max(self.son))
            # nbor
            nbor = f.read_ints()
            self.nbor = nbor.reshape((6,self.noct))
            if not(Silent): print("INFO nbor     ",np.shape(self.nbor),np.min(self.nbor),np.max(self.nbor))
            # octlevel
            self.octlevel = f.read_ints()
            if not(Silent): print("INFO octlevel ",np.shape(self.octlevel),np.min(self.octlevel),np.max(self.octlevel))
            # xoct
            xoct = f.read_reals('d')
            self.xoct = xoct.reshape((3,self.noct))
            if not(Silent): print("INFO xoct     ",np.shape(self.xoct),np.min(self.xoct),np.max(self.xoct)) ###,np.dtype(xoct[0,0])

            # read type gas (composition dependant...)
            if gasmix is None:
                gasmix = 'HI_D_dust'
            if gasmix == ('HI'):
                # for pure HI composition
                if not(Silent): print("-----> gas")
                # velocity
                v = f.read_reals('d')
                v = v.reshape((self.nleaf,3))
                if not(Silent): print("INFO gas v:",np.shape(v),np.amin(v),np.amax(v))
                # nHI
                nHI = f.read_reals('d')
                if not(Silent): print("INFO gas nHI:",np.shape(nHI),np.amin(nHI),np.amax(nHI))
                # dopwidth
                dopwidth = f.read_reals('d')
                if not(Silent): print("INFO gas dopwidth:",np.shape(dopwidth),np.amin(dopwidth),np.amax(dopwidth))
                # boxsize
                [box_size_cm] = f.read_reals('d')
                if not(Silent): print("boxsize [cm] =",box_size_cm)
                f.close()
                # get leaf positions
                xleaf = self.get_leaf_position(Silent)
                # get leaf level
                leaflevel = self.get_leaf_level(Silent)
                # Re-indexing gas mix arrays
                ileaf = np.where(self.son<0)
                icell = np.abs(self.son[ileaf])
                #print(np.shape(icell), np.amin(icell), np.amax(icell))
                ndust = 0
                self.gas = gas(gasmix, nHI[icell-1], dopwidth[icell-1], v[icell-1,:], ndust, xleaf, leaflevel)
            elif gasmix == ('HI_D_dust' or 'HI_dust'):
                # for HI (D) dust composition
                if not(Silent): print("-----> gas")
                # velocity
                v = f.read_reals('d')
                v = v.reshape((self.nleaf,3))
                if not(Silent): print("INFO gas v:",np.shape(v),np.amin(v),np.amax(v))
                # nHI
                nHI = f.read_reals('d')
                if not(Silent): print("INFO gas nHI:",np.shape(nHI),np.amin(nHI),np.amax(nHI))
                # dopwidth
                dopwidth = f.read_reals('d')
                if not(Silent): print("INFO gas dopwidth:",np.shape(dopwidth),np.amin(dopwidth),np.amax(dopwidth))
                # ndust
                ndust = f.read_reals('d')
                if not(Silent): print("INFO gas ndust:",np.shape(ndust),np.amin(ndust),np.amax(ndust))
                # boxsize
                [box_size_cm] = f.read_reals('d')
                if not(Silent): print("boxsize [cm] =",box_size_cm)
                f.close()
                # get leaf positions
                xleaf = self.get_leaf_position(Silent)
                # get leaf level
                leaflevel = self.get_leaf_level(Silent)
                # Re-indexing gas mix arrays
                ileaf = np.where(self.son<0)
                icell = np.abs(self.son[ileaf])
                #print(np.shape(icell), np.amin(icell), np.amax(icell))
                self.gas = gas(gasmix, nHI[icell-1], dopwidth[icell-1], v[icell-1,:], ndust[icell-1], xleaf, leaflevel)
            elif gasmix == 'SiII_dust':
                # for SiII and dust composition
                if not(Silent): print("-----> gas")
                # velocity
                v = f.read_reals('d')
                v = v.reshape((self.nleaf,3))
                if not(Silent): print("INFO gas v:",np.shape(v),np.amin(v),np.amax(v))
                # nHI
                nSiII = f.read_reals('d')
                if not(Silent): print("INFO gas nSiII:",np.shape(nSiII),np.amin(nSiII),np.amax(nSiII))
                # dopwidth
                dopwidth = f.read_reals('d')
                if not(Silent): print("INFO gas dopwidth:",np.shape(dopwidth),np.amin(dopwidth),np.amax(dopwidth))
                # ndust
                ndust = f.read_reals('d')
                if not(Silent): print("INFO gas ndust:",np.shape(ndust),np.amin(ndust),np.amax(ndust))
                # boxsize
                [box_size_cm] = f.read_reals('d')
                if not(Silent): print("boxsize [cm] =",box_size_cm)
                f.close()
                # get leaf positions
                xleaf = self.get_leaf_position(Silent)
                # get leaf level
                leaflevel = self.get_leaf_level(Silent)                
                # Re-indexing gas mix arrays
                ileaf = np.where(self.son<0)
                icell = np.abs(self.son[ileaf])
                #print(np.shape(icell), np.amin(icell), np.amax(icell))
                self.gas = gas(gasmix, nSiII[icell-1], dopwidth[icell-1], v[icell-1,:], ndust[icell-1], xleaf, leaflevel)
            elif gasmix == 'new_ions':
                # for new ions version...
                if not(Silent): print("-----> gas")
                # velocity
                v = f.read_reals('d')
                v = v.reshape((self.nleaf,3))
                if not(Silent): print("INFO gas v:",np.shape(v),np.amin(v),np.amax(v))
                # density /!\ for one element only !!!!
                nion = f.read_reals('d')
                if not(Silent): print("INFO gas nion:",np.shape(nion),np.amin(nion),np.amax(nion))
                # dopwidth
                dopwidth = f.read_reals('d')
                if not(Silent): print("INFO gas dopwidth:",np.shape(dopwidth),np.amin(dopwidth),np.amax(dopwidth))
                # vturb
                vturb = f.read_reals('d')
                if not(Silent): print("INFO gas vturb:",np.shape(vturb),np.amin(vturb),np.amax(vturb))
                # ndust
                ndust = f.read_reals('d')
                if not(Silent): print("INFO gas ndust:",np.shape(ndust),np.amin(ndust),np.amax(ndust))
                # boxsize
                [box_size_cm] = f.read_reals('d')
                if not(Silent): print("boxsize [cm] =",box_size_cm)
                f.close()
                # get leaf positions
                xleaf = self.get_leaf_position(Silent)
                # get leaf level
                leaflevel = self.get_leaf_level(Silent)
                # Re-indexing gas mix arrays
                ileaf = np.where(self.son<0)
                icell = np.abs(self.son[ileaf])
                #print(np.shape(icell), np.amin(icell), np.amax(icell))
                self.gas = gas(gasmix, nion[icell-1], dopwidth[icell-1], v[icell-1,:], ndust[icell-1], xleaf, leaflevel, vturb)
            elif gasmix == 'dust':
                # for new ions version...
                if not(Silent): print("-----> gas")
                # velocity
                v = f.read_reals('d')
                v = v.reshape((self.nleaf,3))
                if not(Silent): print("INFO gas v:",np.shape(v),np.amin(v),np.amax(v))
                # density /!\ for one element only !!!!
                #nion = f.read_reals('d')
                nion = 0
                #if not(Silent): print("INFO gas nion:",np.shape(nion),np.amin(nion),np.amax(nion))
                # dopwidth
                dopwidth = f.read_reals('d')
                if not(Silent): print("INFO gas dopwidth:",np.shape(dopwidth),np.amin(dopwidth),np.amax(dopwidth))
                # vturb
                vturb = f.read_reals('d')
                if not(Silent): print("INFO gas vturb:",np.shape(vturb),np.amin(vturb),np.amax(vturb))
                # ndust
                ndust = f.read_reals('d')
                if not(Silent): print("INFO gas ndust:",np.shape(ndust),np.amin(ndust),np.amax(ndust))
                # boxsize
                [box_size_cm] = f.read_reals('d')
                if not(Silent): print("boxsize [cm] =",box_size_cm)
                f.close()
                # get leaf positions
                xleaf = self.get_leaf_position(Silent)
                # get leaf level
                leaflevel = self.get_leaf_level(Silent)
                # Re-indexing gas mix arrays
                ileaf = np.where(self.son<0)
                icell = np.abs(self.son[ileaf])
                #print(np.shape(icell), np.amin(icell), np.amax(icell))
                self.gas = gas(gasmix, nion, dopwidth[icell-1], v[icell-1,:], ndust[icell-1], xleaf, leaflevel, vturb)
            else:
                #return IOError("mix not defined",gasmix)
                raise NameError("mix not defined",gasmix)
            if not(Silent): print("-----> reading done.")
            if not(Silent): print()


    def get_leaf_position(self,Silent):    
        """ get the leaf cell positions from oct positions"""

        if not(Silent): print("-----> get xleaf")
        ileaf = np.where(self.son<0)
        #print(ileaf)
        # <<< not necessary...
        #ileaf = np.array(ileaf)
        #ileaf = ileaf.reshape(ileaf.size)
        #print(ileaf.dtype)
        #print(ileaf.shape)
        #print(ileaf)
        # >>>
        #icell = np.abs(self.son[ileaf])
        # by definition in rascas
        #ind  = (icell - ncoarse - 1)/noct + 1
        #ioct = icell - ncoarse - ind*noct
        # here icell is ileaf
        ind  = (ileaf - self.ncoarse - 1)//self.noct + 1
        #print(ind)
        ioct = ileaf - self.ncoarse - (ind-1)*self.noct
        #print("ind info: ",np.shape(ind),np.min(ind),np.max(ind))#,np.dtype(ind[0])
        #print("ioct info:",np.shape(ioct),np.min(ioct),np.max(ioct))
        xleaf = np.zeros((3,self.nleaf))
        #print(np.shape(self.octlevel))
        #print(ioct)
        cell_level = self.octlevel[ioct]
        #cell_level = np.empty([self.nleaf],dtype=int)
        #print("cell_level info:",np.shape(cell_level),np.min(cell_level),np.max(cell_level))
        #ioct=ioct.reshape(ioct.size)
        #for i in range(self.nleaf):
        ##    print(i, ioct[i])
        #   cell_level[i] = self.octlevel[ioct[i]-1]
        #print("cell_level info:",np.shape(cell_level),np.min(cell_level),np.max(cell_level))

        # this is the convention used in RASCAS
        offset = np.array((-.5,-.5,-.5, +.5,-.5,-.5, -.5,+.5,-.5, +.5,+.5,-.5, -.5,-.5,+.5, +.5,-.5,+.5, -.5,+.5,+.5, +.5,+.5,+.5),dtype=float).reshape((8,3))
        #print(offset[0,:])
        #print(offset[1,:])
        #print(offset[2,:])
        #print(offset[3,:])
        #print(offset[4,:])
        #print(offset[5,:])
        #print(offset[6,:])
        #print(offset[7,:])

        # cell size
        dx = 0.5**cell_level
        dx = np.reshape(dx,dx.size)
        if not(Silent): print("INFO dx:   ",np.shape(dx), np.amin(dx),np.amax(dx))

        xleaf[0,:] = self.xoct[0,ioct] + offset[ind-1,0]*dx
        xleaf[1,:] = self.xoct[1,ioct] + offset[ind-1,1]*dx
        xleaf[2,:] = self.xoct[2,ioct] + offset[ind-1,2]*dx
        if not(Silent): print("INFO xleaf:",np.shape(xleaf),np.min(xleaf),np.max(xleaf))
    
        return xleaf


    def get_leaf_level(self,Silent):
        """ get the level of leaf cells"""

        if not(Silent): print("-----> get leaf level")
        ileaf = np.where(self.son<0)
        #icell = np.abs(self.son[ileaf])
        ind  = (ileaf - self.ncoarse - 1)//self.noct + 1
        ioct = ileaf - self.ncoarse - (ind-1)*self.noct
        #print("INFO ind:       ",np.shape(ind),np.amin(ind),np.amax(ind))#,np.dtype(ind[0])
        #print("INFO ioct:      ",np.shape(ioct),np.amin(ioct),np.amax(ioct))
        cell_level = np.copy(self.octlevel[ioct])
        cell_level = np.reshape(cell_level,cell_level.size)
        if not(Silent): print("INFO cell_level:",np.shape(cell_level), np.amin(cell_level),np.amax(cell_level))
        return cell_level

