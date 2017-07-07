#!/usr/bin/python

import matplotlib
from matplotlib.pylab import *
import pyfits as pf
from matplotlib.ticker import *
import numpy as np
import matplotlib.font_manager
import matplotlib.cm as cm
from matplotlib import colors
from matplotlib.patches import Ellipse as ellipse
from scipy.optimize import *
from scipy.ndimage.interpolation import *
import seaborn as sns

# Constants and conversion factors
#---------------------------------

# Solar mass (g)
msol = 1.989e+33

# Gravitational constant (cm3 g-1 cm-2)
grav = 6.67259e-8

# Conversion factor: AU to cm
AU2cm = 1.49597871e+13

# Conversion factor: cm to km
cm2km = 1.0e-5

# Speed of light (km)
ckm = 2.99792458e+05

# Disk and star inputs
#---------------------------------

# Distance to source (parsec)
dist = 96.9

# Disk inclination (degrees --> radians)
inc = 36.0*np.pi/180.0

# Position angle (degrees --> radians)
dpa = (144.0 - 90.0)*np.pi/180.0

# Radius r0 (AU) for flaring relation: z = (h0/r0)**psi
r0 = 1.0

# Scale height at r0
h0 = 0

# Flaring index
psi = 1.00

# Mass of central star (g)
mstar = 2.5*msol

# Source velocity (km/s)
source = 5.775

#---------------------------------
# Compute first moment map
#---------------------------------

# Name of FITS image containing first moment map
fits_mom1 = 'CO_346GHz_line.shifted.robust.5.3sig_mom1.fixed.fits'

# Name of FITS image containing the channel maps
fits_chan = 'CO_346GHz_line.shifted.robust.5.fits'

# Read the header from the first moment map
head = pf.getheader(fits_mom1)

# Generate x and y axes: offset position in centimeters
nx    = head['NAXIS1']
xpix  = head['CRPIX1']
xdelt = head['CDELT1']

ny    = head['NAXIS2']
ypix  = head['CRPIX2']
ydelt = head['CDELT2']

# Convert from degrees --> arsec --> AU --> cm
xorig = ((arange(nx) - xpix + 1) * xdelt)*3600
yorig = ((arange(ny) - ypix + 1) * ydelt)*3600

# Source is contained within 4" x 4" - clip the x and y axes extent

xval = xorig[(xorig>=-4.0) & (xorig<=4.0)]*dist*AU2cm
yval = yorig[(yorig>=-4.0) & (yorig<=4.0)]*dist*AU2cm

# Make 2D array containing all (x,y) axis coordinates
xdisk,ydisk = np.meshgrid(xval,yval)


# Read the header from the channel map
head = pf.getheader(fits_chan)

# Create channel axis (Hz)
nv    = head['NAXIS3']
vpix  = head['CRPIX3']
vdelt = head['CDELT3']
vval  = head['CRVAL3']

vel = ((arange(nv) - vpix + 1) * vdelt) + vval

# Extract rest frequency (Hz)
freq  = head['RESTFRQ']

# Convert from frequency (Hz) to LRSK velocity (km/s) using
# rest frequency and source velocity

for i in range(nv):
	vel[i] = (ckm*((freq-vel[i])/freq))-source


# Rx(theta) = | cos(theta) -sin(theta) |
#             | sin(theta)  cos(theta) |

# (x',y') = Rx(theta)*(x,y)
# Note that z-coordinates remain unchanged

# Positive postion angle means clockwise rotation

x0 = np.cos(-dpa)*xdisk - np.sin(-dpa)*ydisk
y0 = np.sin(-dpa)*xdisk + np.cos(-dpa)*ydisk

xdisk = x0
ydisk = y0

# Create a sin(theta) and cos(theta) array where theta is measured from the x-axis

costheta = np.zeros((len(xval),len(yval)),dtype = float)
sintheta = np.zeros((len(xval),len(yval)),dtype = float)

for j in range(len(xval)):
	for i in range(len(yval)):
        	if xval[j] == 0.0 and yval[i] == 0.0:
                	sintheta[i][j] = 0.0
                	costheta[i][j] = 0.0
                else:
                	sintheta[i][j] = xdisk[i][j]/np.sqrt(xdisk[i][j]**2.0 + ydisk[i][j]**2.0)
                	costheta[i][j] = ydisk[i][j]/np.sqrt(xdisk[i][j]**2.0 + ydisk[i][j]**2.0)


#--------------------------------------------------------------------
# Step 2: calculate the z-coordinates for an inclined conical suface
#--------------------------------------------------------------------
# Conical surfaces are simply a series of circles swept out and
# offset and parallel to the plane of symmetry (a geometrically flat disk)

zflat = ydisk*np.tan(inc)

rflat = np.sqrt(xdisk**2.0+ydisk**2.0+zflat**2.0)

# Analytical solution for psi = 1.0
dtop = rflat*(h0/r0)/(1-(h0/r0)*np.tan(inc)*costheta)
dbot = rflat*(h0/r0)/(1+(h0/r0)*np.tan(inc)*costheta)

# Flatten the 2-D arrays
# rinput   = rflat.flatten()/AU2cm
# cosinput = costheta.flatten()
# dtguess  = dtop.flatten()/AU2cm
# dbguess  = dbot.flatten()/AU2cm
# zguess   = zflat.flatten()/AU2cm

# For psi!=1.0, the equation for zprime is a non-linear algebraic equation
#def ftop(z):
#	result = z - (h0/r0)*(rinput+z*np.tan(inc)*cosinput)**psi
#        return result

#def fbot(z):
#	result = z - (h0/r0)*(rinput-z*np.tan(inc)*cosinput)**psi
#        return result

# For psi != 1.0, define the analytical jacobian: diagonals only populated?
#def jtop(z):
#        result = np.zeros((len(z),len(z)))
#        for i in range(len(z)):
#        	result[i][i] =  1 - (h0/r0)*np.tan(inc)*cosinput[i]*psi*(rinput[i]+z[i]*np.tan(inc)*cosinput[i])**(psi-1)
#        return result

#def jbot(z):
#        result = np.zeros((len(z),len(z)))
#        for i in range(len(z)):
#        	result[i][i] =  1 + (h0/r0)*np.tan(inc)*cosinput[i]*psi*(rinput[i]-z[i]*np.tan(inc)*cosinput[i])**(psi-1)
#	return result

# Use scipy.optimize.fsolve to find root of the equation (calls MINPACK's HYBR routines):
#dtop = fsolve(ftop,dtguess,xtol=5.0e-3)
#dbot = fsolve(fbot,dbguess,xtol=5.0e-3)

# Reshape the arrays
dtop = np.reshape(dtop, (len(xval),len(yval)))
dbot = np.reshape(dbot, (len(xval),len(yval)))

# Final scaling by cos(inc)
# dtop = dtop*AU2cm/np.cos(inc)
# dbot = dbot*AU2cm/np.cos(inc)
dtop = dtop/np.cos(inc)
dbot = dbot/np.cos(inc)

# Create the top and bottom cones
ztop = zflat + dtop
zbot = zflat - dbot

#----------------------------------------------
# Step 2: Rosenfeld et al. prescription
#----------------------------------------------

# t^2 * [cos(2i) + cos(2p)] -
#	2 sin^2(p)*[x'^2 sec^2(i) + y'^2 + 2 t x' tan(i)] = 0
#
# a = cos(2i) + cos(2p)
# b = 4 x' sin^2(p) tan(i)
# c = 2 sin^2(p)*[x'^2 sec^2(i) + y'^2]
#
# t = (-b +/- sqrt(b^2 - 4ac))/(2a)

# Angle of flared surface

#phi = np.arctan(h0/r0)

# Quadratic equation coefficients

#a_quad = np.cos(2.0*inc) + np.cos(2.0*phi)
#b_quad = -4.0*ydisk*(np.sin(phi)**2.0)*np.tan(inc)
#c_quad = -2.0*(np.sin(phi)**2.0)*((ydisk**2.0)*(np.cos(inc)**-2.0)+(xdisk**2.0))

# Upper and lower coordinates of the conical surface

#ttop = (-b_quad + np.sqrt(b_quad**2.0 - 4.0*a_quad*c_quad))/(2*a_quad)
#tbot = (-b_quad - np.sqrt(b_quad**2.0 - 4.0*a_quad*c_quad))/(2*a_quad)

#ytop = ydisk/np.cos(inc) + ttop*np.sin(inc)
#xtop = xdisk
#ztop = ttop*np.cos(inc)

#ybot = ydisk/np.cos(inc) + tbot*np.sin(inc)
#xbot = xdisk
#zbot = tbot*np.cos(inc)

#-----------------------------------------------------------
# Step 3: calculate true disk radii projected on x-y plane
#-----------------------------------------------------------

# My prescription uses spherical radius
rtop  = np.sqrt(xdisk**2.0 + ydisk**2.0 + ztop**2.0)
rbot  = np.sqrt(xdisk**2.0 + ydisk**2.0 + zbot**2.0)

# Rosenfeld prescription uses cylindrical radius
#rtop  = np.sqrt(xtop**2.0 + ytop**2.0)
#rbot  = np.sqrt(xbot**2.0 + ybot**2.0)

# Both approaches agree in the limit of small opening angle

#------------------------------------------------
# Step 4: calculate projected keplerian velocity
#------------------------------------------------

vtop  = np.zeros((len(xval),len(yval)),dtype=float)
vbot  = np.zeros((len(xval),len(yval)),dtype=float)
vflat = np.zeros((len(xval),len(yval)),dtype=float)

for j in range(len(xval)):
    for i in range(len(yval)):

               if rtop[i][j]!=0.0:
		vtop[i][j] = cm2km*np.sqrt(grav*mstar/rtop[i][j])*np.sin(inc)*sintheta[i][j]

               if rbot[i][j]!=0.0:
		vbot[i][j] = cm2km*np.sqrt(grav*mstar/rbot[i][j])*np.sin(inc)*sintheta[i][j]

               if rflat[i][j]!=0.0:
               	vflat[i][j] = cm2km*np.sqrt(grav*mstar/rflat[i][j])*np.sin(inc)*sintheta[i][j]
#
#-----------------------------------------------------------------------
# Final step: make a contour map of upper and lower projected velocities
#-----------------------------------------------------------------------

# Reset x and y axes to offset in arcseconds
xval = xval/(dist*AU2cm)
yval = yval/(dist*AU2cm)

# Creat dictionaries containing axes and images
upper_vel = {'image':vtop, 'ra_offset': xval, 'dec_offset': yval}
lower_vel = {'image':vbot, 'ra_offset': xval, 'dec_offset': yval}
flat_vel  = {'image':vflat, 'ra_offset': xval, 'dec_offset': yval}

# Use original x and y axes coordinates for the channel images
img_chan  = {'image':squeeze(pf.getdata(fits_chan)), 'ra_offset':xorig, 'dec_offset':yorig}

# Set spacing between axes labels and tick direction
rcParams['ytick.major.pad'] = 12
rcParams['xtick.major.pad'] = 12

rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

rcParams['xtick.major.size'] = 5
rcParams['xtick.minor.size'] = 2
rcParams['ytick.major.size'] = 5
rcParams['ytick.minor.size'] = 2

rcParams['ytick.labelsize'] = 14
rcParams['xtick.labelsize'] = 14

# Set figure size and create image canvas (in cm)
fig = figure(figsize=[16,8])

# Set seaborn plot styles
sns.set_style("white")
sns.set_style("ticks")
# sns.set_style({"xtick.direction": "in","ytick.direction": "in", "xtick.color" : "white", 'ytick.color' : 'white'})

# Set axes limits
xmin = -4.0
xmax =  4.0
ymin = -4.0
ymax =  4.0

# Choose negative and positive velocity contours to plot
negcontours = ([-0.30,-0.45,-0.60,-0.75,-0.95,-1.05,-1.20,-1.50])
poscontours = ([0.30,0.45,0.60,0.75,0.95,1.05,1.20,1.50])

# Select indices of corresponding channel maps

negchan = np.zeros(len(negcontours))
poschan = np.zeros(len(poscontours))

for i in range(nv-1):
	for j in range(len(negchan)):
		if vel[i] < negcontours[j] and vel[i+1] > negcontours[j]: negchan[j] = i
		if vel[i] < poscontours[j] and vel[i+1] > poscontours[j]: poschan[j] = i

# Set physical range of colour map
cxmin = img_chan['ra_offset'][-1]
cxmax = img_chan['ra_offset'][0]
cymin = img_chan['dec_offset'][0]
cymax = img_chan['dec_offset'][-1]

# Set limits and tics of colorbar
cbmin = -0.1
cbmax = 2.50001
cbtmj = 0.50
cbtmn = 0.10
cbnum = int((cbmax-cbmin)/cbtmn)

# Set colorpalette
#cpal = colors.ListedColormap(sns.cubehelix_palette(cbnum,start=0.5,rot=-0.8,light=0.05,dark=0.95,hue=0.75,reverse=True))
cpal = cm.gist_heat

# Adjust spacing between subplots
subplots_adjust(wspace=0,hspace=0)

# Limit contour extent: within [0:67]
c1 = 12
c2 = 55


for i in range(len(negchan)):

    ax  = fig.add_subplot(2,len(negchan)/2,i+1)
    ax.grid(False)

    ax.set_xlim(xmax,xmin)
    ax.set_ylim(ymin,ymax)

    majorLocator = MultipleLocator(2)
    ax.xaxis.set_major_locator(majorLocator)
    ax.yaxis.set_major_locator(majorLocator)

    minorLocator   = MultipleLocator(0.5)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_minor_locator(minorLocator)

    if i !=len(negchan)/2:
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    else:
        ax.set_xlabel('Relative Right Ascension (arcsec)',x = 0.5, fontsize=14)
        ax.set_ylabel('Relative Declination (arcsec)',y = 0.425, fontsize=14)

        beam_ellipse_color = 'w'
        bmin = head['bmin']*3600.
        bmaj = head['bmaj']*3600.
        bpa  = head['bpa']
        el   = ellipse(xy=[3,-3], width=bmin, height=bmaj, angle=-bpa, \
                fc='k', ec=beam_ellipse_color, fill=True, zorder=10)
        ax.add_artist(el)

    velocity = '+/-%4.2f km/s' % (vel[poschan[i]])
    text(3.5,3.0, velocity, fontsize=14, color = 'w')

    plot([0.0],[0.0], '+', markersize=8, markeredgewidth=1.5, color='w')
    plot([3.5, -3.5],[2.429, -2.5429], ':', color='w')

    # Create rotated images so that the disk major axis runs North to South (i.e. positive and negative velocity components lie on opposite sides of the y-axis).
    negchan_rotated = rotate(img_chan['image'][negchan[i]], 146, reshape = False)
    poschan_rotated = rotate(img_chan['image'][poschan[i]], 146, reshape = False)
#Set the upper half of the negative channel image and the lower half of the positive channel image equal to zero. Add the two images so that contributions from the positive velocity channels are only shown in the upper half of the image, and vice versa. Rotate back to original position angle.
    negchan_rotated[:][:128] = 0
    poschan_rotated[:][128:] = 0
    mirror = rotate(negchan_rotated + poschan_rotated, -146)

    cmap = imshow(mirror,
        	extent=[cxmax,cxmin,cymax,cymin],
            vmin = cbmin,
        	vmax = cbmax,
        	interpolation='bilinear',
		cmap=cpal)

# Generate continuum contour levels and plot contour maps

    # contours = ([negcontours[i]])
    # C1 = contour(-upper_vel['ra_offset'][c1:c2],upper_vel['dec_offset'][c1:c2],
	         # upper_vel['image'][c1:c2,c1:c2],
    #     	 extent=[cxmax,cxmin,cymax,cymin],
	         # levels=contours,
	         # colors='gray',
	         # linewidth=2.0)

    # contours = ([poscontours[i]])
    # C1 = contour(-upper_vel['ra_offset'][c1:c2],upper_vel['dec_offset'][c1:c2],
	         # upper_vel['image'][c1:c2,c1:c2],
    #     	 extent=[cxmax,cxmin,cymax,cymin],
	         # levels=contours,
	         # colors='gray',
	         # linewidth=2.0)

    # contours = ([negcontours[i]])
    # C2 = contour(-lower_vel['ra_offset'][c1:c2],upper_vel['dec_offset'][c1:c2],
	         # lower_vel['image'][c1:c2,c1:c2],
    #     	 extent=[cxmin,cxmax,cymax,cymin],
	         # levels=contours,
	         # colors='white',
	         # linewidth=2.0)

    # contours = ([poscontours[i]])
    # C2 = contour(-lower_vel['ra_offset'][c1:c2],upper_vel['dec_offset'][c1:c2],
	         # lower_vel['image'][c1:c2,c1:c2],
    #     	 extent=[cxmax,cxmin,cymax,cymin],
	         # levels=contours,
	         # colors='white',
	         # linewidth=2.0)

    contours = ([negcontours[i]])
    C2 = contour(-flat_vel['ra_offset'][8:59],flat_vel['dec_offset'][8:59],
	         flat_vel['image'][8:59,8:59],
        	 extent=[cxmax,cxmin,cymax,cymin],
	         levels=contours,
	         colors='black',
	         linewidth=2.0)

    contours = ([poscontours[i]])
    C2 = contour(-flat_vel['ra_offset'][8:59],flat_vel['dec_offset'][8:59],
	         flat_vel['image'][8:59,8:59],
        	 extent=[cxmax,cxmin,cymax,cymin],
	         levels=contours,
	         colors='black',
	         linewidth=2.0)

#	angle=np.array([30,35,40])
#        angle= angle*np.pi/180.0
#
#        y0 = 2.5
#	x0 = 0.5*(y0*np.tan(angle))*cos(inc)
#
#        x1 = np.cos(np.pi/2-dpa)*x0 - np.sin(np.pi/2-dpa)*y0
#        y1 = np.sin(np.pi/2-dpa)*x0 + np.cos(np.pi/2-dpa)*y0
#
#        x2 = np.cos(3*np.pi/2-dpa)*-x0 - np.sin(3*np.pi/2-dpa)*y0
#        y2 = np.sin(3*np.pi/2-dpa)*-x0 + np.cos(3*np.pi/2-dpa)*y0
#
#        for j in range(len(angle)):
#
#                if i < 4:
#                	plot ([0,x1[j]],[0,y1[j]],  ':', linewidth=1, color='k')
#			plot ([0,x2[j]],[0,y2[j]],  ':', linewidth=1, color='k')
#
#                if i >= 4:
#                	plot ([0,x1[j]],[0,y1[j]],  ':', linewidth=1, color='white')
#			plot ([0,x2[j]],[0,y2[j]],  ':', linewidth=1, color='white')

# Create the colorbar
cax = fig.add_axes([0.9, 0.108, 0.02, 0.785])
cbar = colorbar(cmap,cax=cax)
cbar.set_label("Jy/beam",labelpad = 15,fontsize=14)
cbar.set_ticks(np.arange(0,cbmax,cbtmj))
minorLocator   = LinearLocator(cbnum+1)
cbar.ax.yaxis.set_minor_locator(minorLocator)
cbar.update_ticks()


# Save the figure to pdf/eps

fig.savefig('CO_346GHz_line.shifted.robust.5.mirror_channels_flat.png')
show(block=True)