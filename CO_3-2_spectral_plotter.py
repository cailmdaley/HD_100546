import numpy as np
from astropy.io import fits
import matplotlib
from matplotlib.pylab import *
import matplotlib.pyplot as plt
import pyfits as pf
from matplotlib.ticker import *
import matplotlib.font_manager
#import seaborn as sns

#Close previous figure, if open:
#close()

# Name of FITS spectral profile:
fits_spec = 'CO_346GHz_line.shifted.robust.5.spectral_profile.fits'

# Source velocity (km/s)
source = 5.675

# Set conversion factors
cm2km = 1e-5

# Set constants and parameters (in cgs/astronomy units)
ckm  = 2.99792458e5   # cm/km

# Read the header from the observed FITS image
head = pf.getheader(fits_spec)

# Create channel axis (Hz)
nv  = head['NAXIS1']
vdelt = head['CDELT1']
vval  = head['CRVAL1']

vel = ((arange(nv)) * vdelt) + vval

# Extract rest frequency (Hz)
freq  = head['RESTFRQ']

# Convert from frequency (Hz) to LRSK velocity (km/s) using
# rest frequency and source velocity
for i in range(nv):
	vel[i] = (ckm*((freq-vel[i])/freq))-source

# Read in flux densities:
flux_den = fits.getdata(fits_spec)

# Calculate rms of line-free channels:
line_free_chans = np.concatenate((flux_den[:25], flux_den[137:]))
rms = np.sqrt(np.sum(line_free_chans**2)/len(line_free_chans))


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

# Set figure size and create image canvas (in cm)
fig = figure(figsize=[10,10])
ax  = fig.add_subplot(111)

# Set axes limits
xmin = -6.0
xmax =  6.0
ymin =  0.0
ymax =  35.0

ax.set_xlim(xmax,xmin)
ax.set_ylim(ymin,ymax)
ax.grid(False)

# Set x and y labels
ax.set_xlabel('Velocity Relative to Source (km/s)', labelpad=15,fontsize=18)
ax.set_ylabel('Flux Density (Jy)', labelpad=15,fontsize=18)

# Set x and y major and minor tics
x_majorLocator = MultipleLocator(2)
y_majorLocator = MultipleLocator(5)
ax.xaxis.set_major_locator(x_majorLocator)
ax.yaxis.set_major_locator(y_majorLocator)

x_minorLocator   = MultipleLocator(0.25)
y_minorLocator   = MultipleLocator(0.5)
ax.xaxis.set_minor_locator(x_minorLocator)
ax.yaxis.set_minor_locator(y_minorLocator)



#Plot spectral profile and inverse x-axis so that velocity increases left-to-right:
step(vel, flux_den, color='Teal', linestyle='dashed', linewidth=2, where='post')
plot(-vel, flux_den, color = 'k', alpha=0.3, linestyle='dashed', drawstyle='steps-post')
ax.invert_xaxis()

#Add slice along velocity=0:
plot ([0,0],[0,35], ':', color='k')



fig.savefig('CO_346GHz_line.shifted.robust.5.spectral_profile.png')
plt.show()
