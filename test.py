import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
#import pyfits
#from astropy.io import fits 
#import astropy

#NSIDE = 32
#m = np.arange(hp.nside2npix(NSIDE))
#hp.mollview(m, title="Mollview image RING")

map_I = hp.read_map('../PlanckData/COM_CMB_IQU-nilc_1024_R2.02_full.fits',field=0)

#hp.fitsfunc.read_map('./PlanckData/COM_CMB_IQU-nilc_1024_R2.02_full.fits',field=0,h=True)

nside = 1024

#-----------------------------------#
#plot the fullfilled circle

#help(hp.pixelfunc.ang2vec)
#bh: define the center of the disc
vec=hp.pixelfunc.ang2vec(phi=0.,theta=np.pi/2)

#bh: define the radial geodesic arc length
radius = 1.0

#help(hp.query_disc)
disc_pix = hp.query_disc(nside,vec,radius,inclusive=False,fact=4,nest=False)

#bh: erase the cmb info inside the disc
map_I[disc_pix]=1.0e-4

hp.mollview(map_I)

plt.show()

#-----------------------------------#
#plot the unfilled circle

#bh: define a small circle
radius2 = 1.0-0.02

disc_pix2 = hp.query_disc(nside,vec,radius2,inclusive=False,fact=4,nest=False)

#bh: subtract the outer disc by the inner disc, get the ring
ring_pix = np.setdiff1d(disc_pix,disc_pix2)

map_I = hp.read_map('../PlanckData/COM_CMB_IQU-nilc_1024_R2.02_full.fits',field=0)

#bh: erase the cmb info along the ring
map_I[ring_pix]=-3.0e-4

hp.mollview(map_I)

plt.show()

#-----------------------------------#
#plot the fullfilled rectangle

#bh: the definition order of the 4 vertices are important! have to be clock-wise or anticlock-wise!!!
vert1 = hp.pixelfunc.ang2vec(phi=0.,theta=np.pi/3.)
vert2 = hp.pixelfunc.ang2vec(phi=np.pi/3.,theta=np.pi/3.)
vert3 = hp.pixelfunc.ang2vec(phi=np.pi/3.,theta=np.pi/2.)
vert4 = hp.pixelfunc.ang2vec(phi=0.,theta=np.pi/2.)

vertice_coord = [[vert1[0],vert1[1],vert1[2]],[vert2[0],vert2[1],vert2[2]],[vert3[0],vert3[1],vert3[2]],[vert4[0],vert4[1],vert4[2]]]

rect_pix = hp.query_polygon(nside,vertice_coord,inclusive=False,fact=4,nest=False)

map_I = hp.read_map('../PlanckData/COM_CMB_IQU-nilc_1024_R2.02_full.fits',field=0)

map_I[rect_pix]=-2.0e-4

hp.mollview(map_I)

plt.show()

#-----------------------------------#
#plot the unfilled rectangle


#bh: define the inner rectangle
vert_s1 = hp.pixelfunc.ang2vec(phi=0.+0.02,theta=np.pi/3.+0.01)
vert_s2 = hp.pixelfunc.ang2vec(phi=np.pi/3.-0.02,theta=np.pi/3.+0.02)
vert_s3 = hp.pixelfunc.ang2vec(phi=np.pi/3.-0.02,theta=np.pi/2.-0.02)
vert_s4 = hp.pixelfunc.ang2vec(phi=0.+0.02,theta=np.pi/2.-0.02)

vertice_s_coord = [[vert_s1[0],vert_s1[1],vert_s1[2]],[vert_s2[0],vert_s2[1],vert_s2[2]],[vert_s3[0],vert_s3[1],vert_s3[2]],[vert_s4[0],vert_s4[1],vert_s4[2]]]

rect2_pix = hp.query_polygon(nside,vertice_s_coord,inclusive=False,fact=4,nest=False)

#bh: subtract the outer disc by the inner disc, get the ring
rect_ring_pix = np.setdiff1d(rect_pix,rect2_pix)

map_I = hp.read_map('../PlanckData/COM_CMB_IQU-nilc_1024_R2.02_full.fits',field=0)

map_I[rect_ring_pix]=-3.0e-4

hp.mollview(map_I)

plt.show()

quit()

#-----------------------------------#
#gnomonicProj

'''
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


'''
hdulist = pyfits.open('./PlanckData/COM_CMB_IQU-nilc_1024_R2.02_full.fits')

hdulist.info()

print hdulist[0].header
print hdulist[1].header
print hdulist[2].header
'''

#hp.mollview(map_I)

#plt.show()







