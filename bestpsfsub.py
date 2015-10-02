######################################################
########### Aaron J. Juarez, Oct 01, 2015 ############
######################################################
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pyfits
from scipy.ndimage.interpolation import shift,rotate
font = {'family':'serif', 'size':11}
plt.rc('font', **font)

data=np.genfromtxt('centroids.txt',dtype=None,delimiter='\t',
                   names=['filepath','x0','y0'],unpack=1)
x0=data['x0']; y0=data['y0']; filepath = data['filepath']


def slice2d(matrix, start_row, end_row, start_col, end_col):
    return np.array([row[start_col:end_col] for row in matrix[start_row:end_row]])

def compare(frame, other_star):
    resid = []
    for i in range(len(filepath)):
        if other_star in filepath[i]:
            diff = frame - pyfits.getdata(filepath[i])
            sum_residuals = sum(map(sum, diff**2)) #sum over 2d array
            resid.append(sum_residuals)
    return np.array(resid)


stacked_roxs12, stacked_roxs42=0,0
imaj_r12,imaj_r42=[],[]
radz = np.linspace(4,100,10)
f = open('bestmatch.txt','w')
for i in range(len(filepath)):
    print i
    imaj = pyfits.getdata(filepath[i])
    if 'roxs12' in filepath[i]:
        chisq = compare(imaj,'roxs42B'); base=27
    else: #roxs42B
        chisq = compare(imaj,'roxs12'); base=0
    print chisq; print np.shape(chisq)
    chisqind = np.where(chisq==min(chisq))[0][0]
    f.write(filepath[i]+'\t'+filepath[chisqind+base]+'\n')
    imaj -= pyfits.getdata(filepath[chisqind+base])
    dx = x0[i] - 512
    dy = y0[i] - 512
    imaj = shift(imaj,[-dy,-dx])
    head = pyfits.getheader(filepath[i])
    PA = head['PARANG'] + head['ROTPPOSN'] - head['EL'] - head['INSTANGL']
    imaj = rotate(imaj,-PA)
    imaj = slice2d(imaj, len(imaj)/2. - 512, len(imaj)/2. + 512,
                   len(imaj)/2. - 512, len(imaj)/2. + 512)
    if 'roxs12' in filepath[i]:
        stacked_roxs12 += imaj
        imaj_r12.append(imaj)
    else: #roxs42B
        stacked_roxs42 += imaj
        imaj_r42.append(imaj)
f.close()


def show_imaj(imaj, cmap):
    fig = plt.figure(figsize=(9,9))
    ax = fig.add_subplot(111); ax.set_aspect('equal')
    vmin = np.min(imaj); print 'min:',vmin
    vmax = np.max(imaj); print 'max:',vmax
    ax.imshow(imaj,interpolation='nearest',origin='lower',
              cmap=cmap,vmin=vmin,vmax=vmax)
    plt.tight_layout()


imaj_r12 = np.dstack(imaj_r12); median_roxs12 = np.median(imaj_r12, axis=2)
imaj_r42 = np.dstack(imaj_r42); median_roxs42 = np.median(imaj_r42, axis=2)


cmap = cm.bone_r
cmap = cm.cubehelix
imaj_list = [stacked_roxs12,stacked_roxs42,median_roxs12,median_roxs42]
fname = ['matchstk_roxs12','matchstk_roxs42','matchmed_roxs12','matchmed_roxs42']
for i in range(len(imaj_list)):
    show_imaj(imaj_list[i], cmap)
    hdu = pyfits.PrimaryHDU(imaj_list[i])
    hdu.writeto(fname[i]+'.fits', clobber=True)

plt.show()


