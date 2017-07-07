#!/usr/bin/python

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pylab import *
import pyfits as pf
from matplotlib.ticker import *
import numpy as np
import matplotlib.font_manager
import matplotlib.cm as cm
from matplotlib import colors
from matplotlib.patches import Ellipse as ellipse
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

# Generate x and y axes: offset position in arcsec
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

#------------------------------------------------------------------------------
# Step 1: rotate about z axis such that the disk major axis is aligned north to south
#------------------------------------------------------------------------------
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

#-----------------------------------------------------------
# Step 3: calculate true disk radii projected on x-y plane
#-----------------------------------------------------------
rtop  = np.sqrt(xdisk**2.0 + ydisk**2.0 + ztop**2.0)
rbot  = np.sqrt(xdisk**2.0 + ydisk**2.0 + zbot**2.0)

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


#-----------------------------------------------------------------------
# Plot channel maps and isovelocity contours
#-----------------------------------------------------------------------

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

# Use original x and y axes coordinates for the channel images
img_chan = {'image':squeeze(pf.getdata(fits_chan)), 'ra_offset':xorig, 'dec_offset':yorig}

# Set font size and type
#font = {'family' : 'sans-serif',
#	'sans-serif' : ['Helvetica'],
#        'size'   : 14}
#matplotlib.rc('font', **font)

# Set spacing between axes labels and tick direction
rcParams['ytick.major.pad'] = 6
rcParams['xtick.major.pad'] = 6

rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

rcParams['xtick.major.size'] = 5
rcParams['xtick.minor.size'] = 2.5
rcParams['ytick.major.size'] = 5
rcParams['ytick.minor.size'] = 2.5

rcParams['ytick.labelsize'] = 12
rcParams['xtick.labelsize'] = 12

# Set seaborn plot styles
sns.set_style("white")
sns.set_style("ticks")

# Set figure size and create image canvas (in cm)
fig = figure(figsize=[10,10])

# Set axes limits
xmin = -4.0
xmax =  4.0
ymin = -4.0
ymax =  4.0

# Set physical range of colour map
cxmin = img_chan['ra_offset'][0]
cxmax = img_chan['ra_offset'][-1]
cymin = img_chan['dec_offset'][-1]
cymax = img_chan['dec_offset'][0]

# Set limits and tics of colorbar
cbmin = -0.0
cbmax = 2.5005
cbtmaj = 0.5
cbtmin = 0.05
cbnum = int((cbmax-cbmin)/cbtmin)

# Set colorpalette
#cpal = colors.ListedColormap(sns.cubehelix_palette(cbnum,start=0.5,rot=-0.8,light=0.05,dark=0.95,hue=0.75,reverse=True,gamma=1.0))
cpal = cm.gist_heat

# Adjust spacing between subplots
subplots_adjust(wspace=0,hspace=0)

# Set contour levels to show continuum absorption
#rms = 19e-3
#contlev = np.array([-10.0,-3.0])
#contlev = contlev*rms


# Limit isovelocity contour extent: within [0:67]
c1 = 12
c2 = 55

# Loop over channels and plot each panel
for i in range(35):

    chan = 61+i
    velocity = '%4.2f' % (vel[chan])

    ax  = fig.add_subplot(5,7,i+1)

    ax.set_xlim(xmax,xmin)
    ax.set_ylim(ymin,ymax)
    ax.grid(False)

    majorLocator = MultipleLocator(2)
    ax.xaxis.set_major_locator(majorLocator)
    ax.yaxis.set_major_locator(majorLocator)

    minorLocator   = MultipleLocator(1)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_minor_locator(minorLocator)

    if i != 28:
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    else:
        ax.set_xlabel('Relative Right Ascension (arcsec)',x = 0.85,fontsize=12)
        ax.set_ylabel('Relative Declination (arcsec)',y = 0.75,fontsize=12)

        beam_ellipse_color = 'w'
        bmin = head['bmin']*3600.
        bmaj = head['bmaj']*3600.
        bpa  = head['bpa']
        el   = ellipse(xy=[3,-3], width=bmin, height=bmaj, angle=-bpa, \
                fc='w', ec=beam_ellipse_color, fill=True, zorder=10)
        ax.add_artist(el)

    text(3.5,3.0, velocity+' km/s',fontsize=12, color = 'w')

    plot([0.0],[0.0], '+', markersize=6, markeredgewidth=1, color='w')

    cmap = imshow(img_chan['image'][:][:][chan],
            extent=[cxmin,cxmax,cymin,cymax],
            vmin = cbmin,
            vmax = cbmax,
            interpolation='bilinear',
            cmap=cpal)


    #Plot isovelocity contours
    # contours = ([vel[chan]])
    # C1 = contour(-upper_vel['ra_offset'][c1:c2],upper_vel['dec_offset'][c1:c2],
	         # upper_vel['image'][c1:c2,c1:c2],
    #     	 extent=[cxmax,cxmin,cymax,cymin],
	         # levels=contours,
	         # colors='gray',
	         # linewidth=2.0)

    # contours = ([vel[chan]])
    # C1 = contour(-upper_vel['ra_offset'][c1:c2],upper_vel['dec_offset'][c1:c2],
	         # upper_vel['image'][c1:c2,c1:c2],
    #     	 extent=[cxmax,cxmin,cymax,cymin],
	         # levels=contours,
	         # colors='gray',
	         # linewidth=2.0)

    # contours = ([vel[chan]])
    # C2 = contour(-lower_vel['ra_offset'][c1:c2],upper_vel['dec_offset'][c1:c2],
	         # lower_vel['image'][c1:c2,c1:c2],
    #     	 extent=[cxmin,cxmax,cymax,cymin],
	         # levels=contours,
	         # colors='white',
	         # linewidth=2.0)

    # contours = ([vel[chan]])
    # C2 = contour(-lower_vel['ra_offset'][c1:c2],upper_vel['dec_offset'][c1:c2],
	         # lower_vel['image'][c1:c2,c1:c2],
    #     	 extent=[cxmax,cxmin,cymax,cymin],
	         # levels=contours,
	         # colors='white',
	         # linewidth=2.0)

    contours = ([vel[chan]])
    C2 = contour(-flat_vel['ra_offset'][8:59],flat_vel['dec_offset'][8:59],
	         flat_vel['image'][8:59,8:59],
        	 extent=[cxmax,cxmin,cymax,cymin],
	         levels=contours,
	         colors='black',
	         linewidth=2.0)

    contours = ([vel[chan]])
    C2 = contour(-flat_vel['ra_offset'][8:59],flat_vel['dec_offset'][8:59],
	         flat_vel['image'][8:59,8:59],
        	 extent=[cxmax,cxmin,cymax,cymin],
	         levels=contours,
	         colors='black',
	         linewidth=2.0)


        #if (i >=11) & (i <=23):
        #	C1 = contour(xval,yval,
        #		img_chan['image'][:][:][chan],
        #        	levels=contlev,
        #        	colors='gray',
        #        	linewidth='0.2')

# Create and add a colour bar with label and tics
cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])

cbar = colorbar(cmap,cax=cax)
cbar.set_label('Jy/beam',labelpad=-5,fontsize=12)
cbar.set_ticks(np.arange(0,cbmax,cbtmaj))
minorLocator   = LinearLocator(cbnum+1)
cbar.ax.yaxis.set_minor_locator(minorLocator)
cbar.update_ticks()

# Save the figure to pdf/eps

fig.savefig('CO_346GHz_line.shifted.robust.5.channel_maps.png')
plt.show()
