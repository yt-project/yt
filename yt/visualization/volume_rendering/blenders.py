import numpy as np

def enhance(im, stdval=6.0, just_alpha=True):
    if just_alpha:
        nz = im[im>0.0]
        im[:] = im[:]/(nz.mean()+stdval*np.std(nz))
    else:
        for c in range(3):
            nz = im[:,:,c][im[:,:,c]>0.0]
            im[:,:,c] = im[:,:,c]/(nz.mean()+stdval*np.std(nz))
            del nz
    np.clip(im, 0.0, 1.0, im)

def enhance_rgba(im, stdval=6.0):
    nzc = im[:,:,:3][im[:,:,:3]>0.0]
    cmax = nzc.mean()+stdval*nzc.std()


    nza = im[:,:,3][im[:,:,3]>0.0]
    if len(nza) == 0:
        im[:,:,3]=1.0
        amax = 1.0
    else:
        amax = nza.mean()+stdval*nza.std()

    im.rescale(amax=amax, cmax=cmax, inline=True)
    np.clip(im, 0.0, 1.0, im)

