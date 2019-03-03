import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pyfits
#from astropy.io import fits 
#import astropy

#NSIDE = 32
#m = np.arange(hp.nside2npix(NSIDE))
#hp.mollview(m, title="Mollview image RING")

map_I = hp.read_map('../PlanckData/COM_CMB_IQU-nilc_1024_R2.02_full.fits',field=0)

#hp.fitsfunc.read_map('./PlanckData/COM_CMB_IQU-nilc_1024_R2.02_full.fits',field=0,h=True)

#print map_I

hpg = hp.projector.GnomonicProj()

def myvec2pix(x,y,z):
	#(x,y,z) is the cardisian coord of unity sphere
	return hp.pixelfunc.vec2pix(x=x,y=y,z=z,nside=1024) 

#gmap = hp.projector.GnomonicProj().projmap(map_I,myvec2pix,rot=(10,10,5),coord='E')
gmap = hpg.set_proj_plane_info(xsize=400)
gmap = hpg.projmap(map_I,myvec2pix,rot=(10,100,15),coord='E')

#help(hp.projector.GnomonicProj)

print hpg.get_proj_plane_info()
print hpg.get_center(lonlat=True)
print hpg.get_fov()
print gmap.shape

#help(hp.projector.GnomonicProj())

#print hp.projector.GnomonicProj().get_proj_plane_info()

np.savetxt('map_data.dat', gmap)

#import data
cmb_map_data_path = './map_data.dat'
cmb_map_data = np.loadtxt(cmb_map_data_path)

plt.imshow(gmap)

plt.show()

plt.imshow(cmb_map_data)

plt.show()


'''
hdulist = pyfits.open('./PlanckData/COM_CMB_IQU-nilc_1024_R2.02_full.fits')

hdulist.info()

print hdulist[0].header
print hdulist[1].header
print hdulist[2].header
'''

#hp.mollview(map_I)

#plt.show()







