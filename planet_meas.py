######################################################
########### Aaron J. Juarez, Oct 01, 2015 ############
######################################################
import numpy as np
import pyfits, os, fnmatch
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def binary_mask(start_row, end_row, start_col, end_col):
    mask = np.zeros((1024,1024))
    for i in range(len(mask)):
        if i < start_col or i > end_col: pass
        else:
            for j in range(len(mask[i])):
                if j > start_row and j < end_row:
                    mask[j][i] = 1
    return mask

#adapted from http://www.lpl.arizona.edu/~ianc/python/_modules/phot.html#centroid
def centroid(im, mask):
    w = np.ones(im.shape)
    xx = np.arange(im.shape[1])
    yy = np.arange(im.shape[0])
    x,y = np.meshgrid(xx,yy)
    x0 = (x*im*mask*w).sum()/(im*mask*w).sum()
    y0 = (y*im*mask*w).sum()/(im*mask*w).sum()
    return (x0,y0)

def show_imaj(imaj, cmap, c0, ones_mask):
    fig = plt.figure(figsize=(9,9))
    ax = fig.add_subplot(111); ax.set_aspect('equal')
    vmin = np.min(imaj); print 'min:',vmin
    vmax = np.max(imaj); print 'max:',vmax
    ax.plot(c0[0],c0[1],'b*',c='gold',ms=9,alpha=0.5)
    ax.imshow(imaj,interpolation='nearest',origin='lower',
              cmap=cmap,vmin=vmin,vmax=vmax)
#    ax.imshow(np.array(ones_mask),interpolation='nearest',origin='lower',
#              cmap=cmap,vmin=vmin,vmax=vmax,alpha=0.4)
    plt.tight_layout()


filenom=[]
for file in os.listdir('.'):
    if fnmatch.fnmatch(file, '*med*.fits'):
        filenom.append(file)
#fig 12, 11 (median_roxs42.fits) -- can't see planets
#fig 10 (median_roxs12.fits) -- blurry
filenom = np.delete(filenom,[6,7]) # ignore these files


cmap = cm.bone_r
cmap = cm.cubehelix
f = open('planets.txt','w')
# http://www2.keck.hawaii.edu/inst/nirc2/genspecs.html
pixscale = 0.009952 #arcsec/pixel (+/- 0.00005"); 10x10 arcsec fov (narrow camera)
for i in range(len(filenom)):
    print i, filenom[i]
    imaj = pyfits.getdata(filenom[i])
    if 'roxs12' in filenom[i]:
        ones_mask = binary_mask(670,700,470,500)
        c0 = centroid(imaj, ones_mask)
#        show_imaj(imaj,cmap,c0,ones_mask)
        projsep = np.sqrt((c0[0]-512)**2+(c0[1]-512)**2) * pixscale
        PA = np.arctan2(c0[1]-512,c0[0]-512)*(180/np.pi)+270
        if PA > 360: PA -= 360
        f.write(filenom[i]+'\t'+str(projsep)+'\t'+str(PA)+'\n')
    else: #roxs42
        ones_mask = binary_mask(495,525,615,645)
        c0 = centroid(imaj, ones_mask)
#        show_imaj(imaj,cmap,c0,ones_mask)
        projsep = np.sqrt((c0[0]-512)**2+(c0[1]-512)**2) * pixscale
        PA = np.arctan2(c0[1]-512,c0[0]-512)*(180/np.pi)+270
        if PA > 360: PA -= 360
        f.write(filenom[i]+'-1\t'+str(projsep)+'\t'+str(PA)+'\n')

        ones_mask = binary_mask(455,480,540,565)
        c0 = centroid(imaj, ones_mask)
#        show_imaj(imaj,cmap,c0,ones_mask)
        projsep = np.sqrt((c0[0]-512)**2+(c0[1]-512)**2) * pixscale
        PA = np.arctan2(c0[1]-512,c0[0]-512)*(180/np.pi)+270
        if PA > 360: PA -= 360
        f.write(filenom[i]+'-2\t'+str(projsep)+'\t'+str(PA)+'\n')


f.close()
plt.show()

