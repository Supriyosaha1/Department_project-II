
import numpy as np

def make_map(lmax,positions,levels,features,xmin,xmax,ymin,ymax,zmin,zmax,view_dir='z'):
    """
    From Harley Katz
    Adapted by Leo Michel-Dansac
    (March 2023)
    """
    
    npix = 2**lmax
    pixel_positions = positions*(2**lmax)
    dx = 1./2.**levels
    # define ranges
    # Get the ranges for the image in pixel pos...
    x_l = xmin*npix
    x_h = xmax*npix
    y_l = ymin*npix
    y_h = ymax*npix
    z_l = zmin*npix
    z_h = zmax*npix

    if view_dir == 'x':
        im_range = ((y_l,y_h),(z_l,z_h))
        i1 = 1
        i2 = 2
        i3 = 0
        zz_l = x_l
        zz_h = x_h
    if view_dir == 'y':
        im_range = ((x_l,x_h),(z_l,z_h))
        i1 = 0
        i2 = 2
        i3 = 1
        zz_l = y_l
        zz_h = y_h
    if view_dir == 'z':
        im_range = ((x_l,x_h),(y_l,y_h))
        i1 = 0
        i2 = 1 
        i3 = 2
        zz_l = z_l
        zz_h = z_h

    l_pix_per_level = [int(npix/(2.**(lmax-l))) for l in range(1,lmax+1)]

    image = np.zeros((npix,npix))
    all_images = []
    for i,l in enumerate(range(1,lmax+1)):
        if l_pix_per_level[i] < 1:
            continue

        # select at least one slice of cells at each level
        zmin_l = (0.5-(1./2**l)) * (2**lmax)
        zmax_l = ((1./2**l) + 0.5) * (2**lmax)
        zmin_l = min(zmin_l,zz_l)
        zmax_l = max(zmax_l,zz_h)
        #filt = levels == l
        filt = (levels == l) & (pixel_positions[i3,:] >= zmin_l) & (pixel_positions[i3,:] <= zmax_l)

        H, _, _ = np.histogram2d(pixel_positions[i1,:][filt],
                                 pixel_positions[i2,:][filt],
                                 bins=l_pix_per_level[i],
                                 range=im_range,weights=features[filt]*dx[filt])

        all_images.append(H)
        if l < lmax:
            up_samp = int(2**(lmax-l))
            H = H.repeat(up_samp, axis=1).repeat(up_samp, axis=0)

        image += H

    return image,all_images


def lmax_map(lmax,positions,levels,xmin,xmax,ymin,ymax,zmin,zmax,view_dir='z'):
    """
    From Harley Katz
    Adapted by Leo Michel-Dansac
    (March 2023)
    """
    
    npix = 2**lmax
    #print('npix = ',npix)
    pixel_positions = positions*(2**lmax)

    #print(np.min(pixel_positions[0,:]), np.max(pixel_positions[0,:]))
    #print(np.min(pixel_positions[1,:]), np.max(pixel_positions[1,:]))
    #print(np.min(pixel_positions[2,:]), np.max(pixel_positions[2,:]))
    # define ranges
    # Get the ranges for the image in pixel pos...
    x_l = xmin*npix
    x_h = xmax*npix
    y_l = ymin*npix
    y_h = ymax*npix
    z_l = zmin*npix
    z_h = zmax*npix

    if view_dir == 'x':
        im_range = ((y_l,y_h),(z_l,z_h))
        i1 = 1
        i2 = 2
        i3 = 0
        zz_l = x_l
        zz_h = x_h
    if view_dir == 'y':
        im_range = ((x_l,x_h),(z_l,z_h))
        i1 = 0
        i2 = 2
        i3 = 1
        zz_l = y_l
        zz_h = y_h
    if view_dir == 'z':
        im_range = ((x_l,x_h),(y_l,y_h))
        i1 = 0
        i2 = 1 
        i3 = 2
        zz_l = z_l
        zz_h = z_h
    
    l_pix_per_level = [int(npix/(2.**(lmax-l))) for l in range(1,lmax+1)]

    image = np.zeros((npix,npix))
    all_images = []
    for i,l in enumerate(range(1,lmax+1)):
        if l_pix_per_level[i] < 1:
            continue

        # select at least one slice of cells at each level
        zmin_l = (0.5-(1./2**l)) * (2**lmax)
        zmax_l = ((1./2**l) + 0.5) * (2**lmax)
        zmin_l = min(zmin_l,zz_l)
        zmax_l = max(zmax_l,zz_h)
        #print(zmin_l,zmax_l)
        
        #filt = levels == l
        filt = (levels == l) & (pixel_positions[i3,:] >= zmin_l) & (pixel_positions[i3,:] <= zmax_l)
        #print(np.shape(levels[filt]))
        
        H, _, _ = np.histogram2d(pixel_positions[i1,:][filt],
                                 pixel_positions[i2,:][filt],
                                 bins=l_pix_per_level[i],
                                 range=im_range,density=False)

        H[H>0] = l
        all_images.append(H)
        if l < lmax:
            up_samp = int(2**(lmax-l))
            H = H.repeat(up_samp, axis=1).repeat(up_samp, axis=0)

        image = np.maximum(image,H)

    return image,all_images
