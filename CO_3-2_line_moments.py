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
#import seaborn as sns
from my_colormaps import *

# Name of FITS images
#fits_cont = './cont_sb2_final.fixed.image.fits'
fits_mom0 = 'CO_346GHz_line.shifted.robust.5.uncut_mom0.fits'
fits_mom1 = 'CO_346GHz_line.shifted.robust.5.3sig_mom1.fixed.fits'

# Read the header from the observed FITS continuum image: axes should be
# same for the lines
head = pf.getheader(fits_mom0)

# Generate x and y axes: offset position in arcsec

nx = head['NAXIS1']
xpix = head['CRPIX1']
xdelt = head['CDELT1']

ny = head['NAXIS2']
ypix = head['CRPIX2']
ydelt = head['CDELT2']

# Convert from degrees to arcsecs
xval = ((arange(nx) - xpix + 1) * xdelt) * 3600
yval = ((arange(ny) - ypix + 1) * ydelt) * 3600

# Creat dictionaries containing axes and images
#img_cont = {'image':squeeze(pf.getdata(fits_cont)), 'ra_offset':xval, 'dec_offset':yval}
img_mom0 = {'image': squeeze(pf.getdata(fits_mom0)),
            'ra_offset': xval, 'dec_offset': yval}
img_mom1 = {'image': squeeze(pf.getdata(fits_mom1)),
            'ra_offset': xval, 'dec_offset': yval}

masked_mom1 = np.ma.array(img_mom1['image'], mask=np.isnan(img_mom1['image']))

# Set font size and type
# font = {'family' : 'sans-serif',
#	'sans-serif' : ['Helvetica'],
#        'size'   : 20}
#matplotlib.rc('font', **font)

# Set spacing between axes labels and tick direction
rcParams['ytick.major.pad'] = 12
rcParams['xtick.major.pad'] = 12

rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

rcParams['xtick.major.size'] = 10
rcParams['xtick.minor.size'] = 5
rcParams['ytick.major.size'] = 10
rcParams['ytick.minor.size'] = 5

rcParams['ytick.labelsize'] = 18
rcParams['xtick.labelsize'] = 18

rcParams['lines.linewidth'] = 1

# Set seaborn plot styles
# sns.set_style("white")
# sns.set_style("ticks")

# Set figure size and create image canvas (in cm)
fig = figure(figsize=[10, 10])
ax = fig.add_subplot(111)

# Set axes limits
xmin = -4.0
xmax = 4.0
ymin = -5.0
ymax = 4.0

ax.set_xlim(xmax, xmin)
ax.set_ylim(ymin, ymax)
ax.grid(False)

# Set x and y labels

ax.set_xlabel('Relative Right Ascension (arcsec)', labelpad=15, fontsize=18)
ax.set_ylabel('Relative Declination (arcsec)', labelpad=15, fontsize=18)

# Set x and y major and minor tics
majorLocator = MultipleLocator(1)
ax.xaxis.set_major_locator(majorLocator)
ax.yaxis.set_major_locator(majorLocator)

minorLocator = MultipleLocator(0.1)
ax.xaxis.set_minor_locator(minorLocator)
ax.yaxis.set_minor_locator(minorLocator)

# Set physical range of colour map
cxmin = img_mom0['ra_offset'][0]
cxmax = img_mom0['ra_offset'][-1]
cymin = img_mom0['dec_offset'][-1]
cymax = img_mom0['dec_offset'][0]

# Set limits and tics of colorbar - velocity scale
cbmin = -3.00
cbmax = 3.001
cbtmj = 0.50
cbtmn = 0.10
cbnum = int((cbmax - cbmin) / cbtmn)

# Set colorpalette
#cpal = colors.ListedColormap(sns.color_palette("coolwarm", cbnum, desat=0.7))
cpal = casa_colors()
cpal.set_bad('w', 1.)

# Plot moment 1 map as a colour map
cmap = imshow(masked_mom1,
              extent=[cxmin, cxmax, cymin, cymax],
              vmin=cbmin,
              vmax=cbmax,
              interpolation='bilinear',
              cmap=cpal)

# Create the colorbar
cbar = colorbar(cmap, shrink=0.8)
cbar.set_label("km/s", labelpad=15, fontsize=18)
cbar.set_ticks(np.arange(cbmin, cbmax, cbtmj))
minorLocator = LinearLocator(cbnum + 1)
cbar.ax.yaxis.set_minor_locator(minorLocator)
cbar.update_ticks()

# Set parameters for contour maps
# cont_rms = 6.3e-04                      # continuum rms noise in Jy/beam
line_rms = 0.055                     # line rms noise in Jy/beam km/s

# Rescale images by sigma
#img_cont['image'] = img_cont['image']/cont_rms
img_mom0['image'] = img_mom0['image'] / line_rms

# cont_cont = [3]                         # contour levels for continuum
line_cont = [3, 10, 30, 100, 150]           # contour levels for line
# col_cont = 'k'                          # contour colours for continuum
col_line = 'k'                       # contour colours for line
cont_width = 2.0                        # contour line width

# Generate continuum contour levels and plot contour map
#contlev = array(cont_cont)
# C1 = contour(ra,dec,
#        img_cont['image'],
#        levels=contlev,
#        colors=col_cont,
#        linestyles='dashed')

# Generate line contour levels and plot contour map
contlev = array(line_cont)
C1 = contour(ra, dec,
             img_mom0['image'],
             levels=contlev,
             colors=col_line,
             linewidth=cont_width,
             linestyles='dashed')

# clabel(C1,inline=True,fontsize=12,fmt="%4.0f")

# Overplot the beam ellipse
beam_ellipse_color = 'k'
bmin = head['bmin'] * 3600.
bmaj = head['bmaj'] * 3600.
bpa = head['bpa']
el = ellipse(xy=[3.5, -4.5], width=bmin, height=bmaj, angle=-bpa,
             fc='k', ec=beam_ellipse_color, fill=True, zorder=10)
ax.add_artist(el)

# Plot the scale bar
#plot([-3,-4],[-4.5,-4.5], '-', linewidth=3, color='k')

# Plot a cross at the source position
plot([0.0], [0.0], '+', markersize=12, markeredgewidth=2, color='k')

# Add figure labels

#figtext(0.61, 0.24, '160 AU', size = 18, color = 'k',backgroundcolor = 'w')

# Plot slice along disk major axis
plot([3.5, -3.5], [2.5429, -2.5429], ':', color='k')
plot([-2.5429, 2.9062], [3.5, -4.0], ':', color='k')

# Save the figure to pdf/eps

fig.savefig('CO_346GHz_line.shifted.robust.5.moments.png')
show(block=True)
