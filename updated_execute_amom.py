from matplotlib.patches import Ellipse as ellipse
from my_colormaps import casa_colors
from matplotlib.pylab import *
from amom import *
import matplotlib.cm as cm
import matplotlib.colors as colors
import pyfits as pf
import seaborn as sns
import itertools

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Grid Parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Parameters for all models
mstar    = 2.5
aspects  = np.array([0.03])
cones    = ['lower']
# aspects = np.array([0.0,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5])
# cones   = ['lower','upper']

# Unwarped disk Parameters
dpas     = np.array([144.]) + 90.
incls    = np.array([36.0])
# dpas   = np.arange(130,160,2) + 90
# incls  = np.arange(30,60,2)

# Warped disk Parameters
warp=False
pa0s     = np.array([114.0]) + 90.0
pa_outs  = np.array([144.0]) + 90.0
i0s      = np.array([56.0])
i_outs   = np.array([36.0])
r0s      = np.array([40.0])
grids    = ['inner']
# pa0s      = np.arange(60,131,3) + 90
# pa_outs   = np.array([144.0]) + 90.0
# i0s       = np.arange(36,77,4)
# i_outs    = np.array([36.0])
# r0s       = np.arange(30, 90, 10)

# Blob Parameters
blob=True

phi_sts    = np.array([255.0])
phi_ends   = np.array([285.0])
r_sts      = np.array([0.])
r_ends     = np.array([80.0])
v_blobs  = np.array([-5.0])

# phi_sts  = np.arange(240, 281, 5)
# phi_ends = np.arange(280,  311, 5)
# r_sts    = np.array([0.])
# # r_ends   = np.array([15, 20, 25,30, 40, 50, 60, 70, 80, 90, 100])
# r_ends  = np.arange(30, 101, 10)
# v_blobs  = np.arange(-14, -2.1, 3)

# Fit Plot x and y axes:
xparam   = phi_sts
yparam   = phi_ends #still have to change np.where("____"==xparam) on lines 142-146


# Read in header from FITS first moment map
head   = pf.getheader('../Plotting/CO_346GHz_line.shifted.robust.5.3sig_mom1.fixed.fits')

# Extract source RA, Dec, and the size of the image in AU from header
ra0 = head['CRVAL1']
dec0 = head['CRVAL2']
dpc    = 96.9
sizeau = abs(head['CDELT1'])*head['NAXIS1']*3600.*dpc

# Generate x and y axes: offset position in arcsec
nx    = head['NAXIS1']
xpix  = head['CRPIX1']
xdelt = head['CDELT1']

ny    = head['NAXIS2']
ypix  = head['CRPIX2']
ydelt = head['CDELT2']

# Convert from degrees to arcsecs
xval = ((arange(nx) - xpix + 1) * xdelt)*3600
yval = ((arange(ny) - ypix + 1) * ydelt)*3600

# Set axes limits in arcseconds
xmin = -4.0
xmax =  4.0
ymin = -4.0
ymax =  4.0

# Initiate variables for interpolation options
N = 'none'
B = 'bilinear'

# Set physical range of fit and colour maps
cxmin = xval[0]
cxmax = xval[-1]
cymin = yval[-1]
cymax = yval[0]

# Set limits and tics of colorbar - velocity scale
cbmin = np.array([-4, -4, -1., -7.0])
cbmax = np.array([ 4.0001,  4.0001,  1.0001,  7.0001])
cbtmj = [1.0, 1.0, 0.25, 1.0]
cbtmn = [0.25, 0.25, 0.1, 0.5]
cbnum = ((cbmax-cbmin)/cbtmn).astype(int)
vel_res_ticks = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]


# Return model and goodness of fit statsitics
def create_model(image, dpc=0., mstar=1., incl=0., dpa=0., aspect=0., cone=0., warp=False, grid='inner', r0=0., pa0=0., pa_out=0., i0=0., i_out=0., blob=False, phi_st=0., phi_end=0., r_st=0., r_end=0., v_blob=0.):
    # Create model first moment map grid
    mom1 = get_mom1(npix=head['NAXIS1'], sizeau=sizeau, dpc=dpc, mstar=mstar, incl=incl, dpa=dpa, aspect=aspect, cone=cone, warp=warp, grid=grid, r0=r0, pa0=pa0, pa_out=pa_out, i0=i0, i_out=i_out, blob=blob, phi_st=phi_st, phi_end=phi_end, r_st=r_st, r_end=r_end, v_blob=v_blob,  x=oim['x'], y=oim['y'], ctrpix=0.)

    # Convolve model with data
    cmom1 = mom1_conv(mom1['mom1'], fwhm=[head['bmaj']*3600., head['bmin']*3600], pa=360.-head['bpa'], pscale=head['CDELT1']*3600.)

    # Create residuals; change the sign of residuals where the data velocity is negative so that residuals can be interpreted consistently across the entire first moment map.
    resid = oim['image'] - cmom1
    resid[np.where(oim['image'] <= 0.0)] *= -1.0

    # Create create residuals that are in units of velocity resolution (.22km/s):
    resid_vel_res = resid/0.22

    # Create dictionary containing axes and images
    images = {'image':[oim['image'], cmom1, resid, resid_vel_res], 'ra_offset':xval, 'dec_offset':yval}


    # Remove nan values so that arithmetic operations can be performed.
    resid_nonans = np.copy(resid)
    nans = np.isnan(resid_nonans)
    resid_nonans[(np.where(nans == True))] = 0

    # Calculate goodness of fit (peak, peak location, and quadratic sum):
    res_peak = np.min(resid_nonans) if np.max(abs(resid_nonans)) == abs(np.min(resid_nonans)) else np.max(resid_nonans)
    peak_pixel_y, peak_pixel_x = np.where(resid_nonans==res_peak)
    peak_offset_y, peak_offset_x = yval[peak_pixel_y][0], xval[peak_pixel_x][0]

    res_quadsum = np.sum(resid_nonans**2)
    res_quadsum_out = np.sum(resid_nonans[np.where( mom1['rdp'] > 150*au)]**2)

    resid_peak_array[np.where(phi_end==yparam), np.where(phi_st==xparam)] = res_peak
    peak_offset_array[np.where(phi_end==yparam), np.where(phi_st==xparam)] = [peak_offset_x, peak_offset_y]

    resid_sum_array[np.where(phi_end==yparam), np.where(phi_st==xparam)] = res_quadsum
    resid_sum_out_array[np.where(phi_end==yparam), np.where(phi_st==xparam)] = res_quadsum_out


    #Print parameters and goodness of fit statistics
    print('~'*61)
    # print('Transition Radius: {:>12}       Grid: {:>30}'.format(r0, grid))
    # print('Outer Inclination: {:>12}       Outer Position Angle: {:>14}'.format(i_out, pa_out-90.0))
    # print('Inner Inclination: {:>12}       Inner Position Angle: {:>14}'.format(i0, pa0-90.0))
    print('Inner Radius: {:>16}       Outer Radius: {:>30}'.format(r_st, r_end))
    print('Start Angle: {:>16}       End Angle: {:>17}'.format(phi_st, phi_end))
    print('Aspect: {:>23}       Cone: {:>30}'.format(aspect, cone))
    print('Blob Velocity: {:>16}'.format(v_blob))
    print
    print
    print('Residual Peak: {0:>16.4f}       Residual Square Sum: {1:>15.4f}'.format(res_peak, res_quadsum))
    return images, [res_peak, peak_offset_x, peak_offset_y, res_quadsum, res_quadsum_out]


def model_plot(images, xlabel, ylabel, x_majmin, y_majmin, extent,interpolation, figsize=[5,5], cpal=casa_colors(), rows=1, cols=1, label_fig=0, xtop=False, cbar_each=False, fig_titles=None, cbar_labels=None, save=None, show=False, fit_stats=None):

    sns.set_style('ticks')
    sns.set_context('talk')

    # Set spacing between axes
    rcParams['ytick.major.pad'] = 5
    rcParams['xtick.major.pad'] = 5

    rcParams['xtick.direction'] = 'in'
    rcParams['ytick.direction'] = 'in'

    rcParams['xtick.major.size'] = 10
    rcParams['xtick.minor.size'] = 5
    rcParams['ytick.major.size'] = 10
    rcParams['ytick.minor.size'] = 5

    rcParams['ytick.labelsize'] = 10
    rcParams['xtick.labelsize'] = 10

    rcParams['lines.linewidth'] = 1

    # Set figure size and create image canvas (in cm)
    fig = figure(figsize=figsize)

    # Set colorpalette
    cpal.set_bad('w',1.)

    # Adjust spacing between subplots
    if cbar_each==False:
        subplots_adjust(wspace=0,hspace=0)


    for i in range((rows-1)*cols):
        ax  = fig.add_subplot(rows,cols,i+1)

        #Set spatial extent of plot
        ax.set_xlim(xmax,xmin)
        ax.set_ylim(ymin,ymax)


        #Set tick marks
        xmajorLocator = MultipleLocator(x_majmin[0])
        xminorLocator = MultipleLocator(x_majmin[1])
        ymajorLocator = MultipleLocator(y_majmin[0])
        yminorLocator = MultipleLocator(y_majmin[1])

        ax.xaxis.set_major_locator(xmajorLocator)
        ax.xaxis.set_minor_locator(xminorLocator)
        ax.yaxis.set_major_locator(ymajorLocator)
        ax.yaxis.set_minor_locator(yminorLocator)

        # Set labels
        if i==label_fig:
            if xtop:
                ax.xaxis.set_label_position('top')
                ax.xaxis.tick_top()
            ax.set_ylabel(ylabel,y = 0.50,fontsize=11)
            ax.set_xlabel(xlabel,fontsize=11)

            ax.yaxis.labelpad=3
            ax.xaxis.labelpad=3

        if cbar_each==True:
            if fig_titles is not None:
                ax.set_title(fig_titles[i], fontsize=13, y=1.02)

        else:
            ax.set_xticklabels([])
            ax.set_yticklabels([])

        # Create beam
        beam_ellipse_color = 'k'
        bmin = head['bmin']*3600.
        bmaj = head['bmaj']*3600.
        bpa  = head['bpa']
        el   = ellipse(xy=[3,-3], width=bmin, height=bmaj, angle=-bpa, \
                fc='0.75', ec=beam_ellipse_color, fill=True, zorder=10)
        ax.add_artist(el)

        #Add cross at source position
        ax.plot([0.0],[0.0], '+', markersize=10, markeredgewidth=2, color='k')


        if i==3:
            norm = colors.BoundaryNorm(vel_res_ticks, cpal.N)
        else: norm=None

        cmap = imshow(images[i],
                extent =extent,
                vmin   = cbmin[i],
                vmax   = cbmax[i],
                norm=norm,
                origin ='upper',
                interpolation=interpolation[i],
                cmap=cpal)

        if cbar_each == True:
            if i!=3:
                cbar = colorbar(cmap)
                cbar.set_label(cbar_labels[i], labelpad=3, fontsize=10)
                cbar.set_ticks(np.arange(cbmin[i],cbmax[i],cbtmj[i]))
                minorLocator   = LinearLocator(cbnum[i]+1)
                cbar.ax.yaxis.set_minor_locator(minorLocator)
                cbar.update_ticks()
            else:
                cbar = colorbar(cmap=cpal, norm=norm, spacing='proportional',ticks=vel_res_ticks, boundaries=vel_res_ticks)
                cbar.set_label(cbar_labels[i], labelpad=3, fontsize=10)

    # Create and add a colour bar with label and ticks
    if cbar_each == False:
        cax = fig.add_axes([0.91, 0.1, 0.03, 0.8])

        cbar = colorbar(cmap,cax=cax)
        cbar.set_label(cbar_labels,labelpad = 1,fontsize=10)
        cbar.set_ticks(np.arange(cbmin,cbmax,cbtmj))
        minorLocator   = LinearLocator(cbnum+1)
        cbar.ax.yaxis.set_minor_locator(minorLocator)
        cbar.update_ticks()

    #Create histogram
    ax  = fig.add_subplot(rows,cols,i+2)

    ax.set_ylabel('Number of Pixels', fontsize=11)
    ax.set_xlabel('Velocity Channel (km/s)', fontsize=11)
    ax.set_title('Residual Pixel Histogram', y=1.02, fontsize=13)

    ax.yaxis.labelpad=3
    ax.xaxis.labelpad=5


    bins = np.arange(-15*0.11, 15*0.11+.001, 0.11)
    hist(images[2][np.invert(np.isnan(images[2]))], bins=bins, log=True)
    ax.axvline(-0.22, ls='--', c='k')
    ax.axvline(0.22, ls='--', c='k')

    plt.tight_layout()


    # Add goodness of fit estimates
    if fit_stats:
        figtext(0.55, 0.325, r'Peak Residual$={0:>10.4f}$'.format(fit_stats[0]), fontsize=15)
        figtext(0.55, 0.31, r'Peak Coordinates$=({0},{1})$'.format(fit_stats[1], fit_stats[2]), fontsize=15)

        figtext(0.55, 0.25, r'Residual Quadratic Sum$= {0:.4f}$'.format(fit_stats[3]), fontsize=15)
        figtext(0.55, 0.235, r'Quadratic Sum Outside 150 AU$= {0:.4f}$'.format(fit_stats[4]), fontsize=15)

    if save:
        savefig(save)
    if show:
        plt.show()


def fit_plot(images, xlabel, ylabel, x_majmin, y_majmin, save=None, show=False, figsize=[5,5], cpal=casa_colors(), rows=1, cols=1, label_fig=0, xtop=False, cbar_label=None, cbar_each=False, vmin=None, vmax=None, bilin=None, fit_stats=None):

    # Set spacing between axes
    rcParams['ytick.major.pad'] = 5
    rcParams['xtick.major.pad'] = 5

    rcParams['xtick.direction'] = 'in'
    rcParams['ytick.direction'] = 'in'

    rcParams['xtick.major.size'] = 10
    rcParams['xtick.minor.size'] = 5
    rcParams['ytick.major.size'] = 10
    rcParams['ytick.minor.size'] = 5

    rcParams['ytick.labelsize'] = 10
    rcParams['xtick.labelsize'] = 10

    rcParams['lines.linewidth'] = 1

    # Set figure size and create image canvas (in cm)
    fig = figure(figsize=figsize)

    # Set colorpalette
    cpal.set_bad('w',1.)


    #Set x,y limits
    # xmin = xparam[0] - 90.
    # xmax = xparam[-1] + (xparam[1]-xparam[0])- 90.
    xmin = xparam[0]
    xmax = xparam[-1] + (xparam[1]-xparam[0])
    ymin = yparam[0]
    ymax = yparam[-1] + (yparam[1]-yparam[0])


    for i in range(rows*cols):
        ax  = fig.add_subplot(rows,cols,i+1)

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.grid(False)

        xmajorLocator = MultipleLocator(x_majmin[0])
        xminorLocator = MultipleLocator(x_majmin[1])
        ymajorLocator = MultipleLocator(y_majmin[0])
        yminorLocator = MultipleLocator(y_majmin[1])

        ax.xaxis.set_major_locator(xmajorLocator)
        ax.xaxis.set_minor_locator(xminorLocator)
        ax.yaxis.set_major_locator(ymajorLocator)
        ax.yaxis.set_minor_locator(yminorLocator)

        xticks(rotation=70)

        if i == label_fig:
            ax.set_ylabel(ylabel,y = 0.50,fontsize=11)
            ax.set_xlabel(xlabel,fontsize=12)


        cmap = imshow(images[i],
                extent=[xmin,xmax,ymin,ymax],
                vmin=vmin[i],
                vmax=vmax[i],
                origin='lower',
                interpolation='none',
                cmap=cpal)

        if cbar_each == True:
            cbar = colorbar(cmap)
            cbar.set_label(cbar_label[i], labelpad=3, fontsize=10)

    # Create and add a colour bar with label and ticks
    if cbar_each == False:
        cax = fig.add_axes([0.91, 0.1, 0.03, 0.8])

        cbar = colorbar(cmap,cax=cax)
        cbar.set_label(cbar_label,labelpad = 1,fontsize=10)
        cbar.set_ticks(np.arange(cbmin,cbmax,cbtmj))
        minorLocator   = LinearLocator(cbnum+1)
        cbar.ax.yaxis.set_minor_locator(minorLocator)
        cbar.update_ticks()

    # Add goodness of fit estimates
    if fit_stats is not None:
        best_fit = np.min(images[0])
        incl_ind, pa_ind = np.where(images[0]==best_fit)
        min_incl, min_pa = yparam[incl_ind][0], xparam[pa_ind][0]
        x_offset, y_offset = fit_stats[incl_ind, pa_ind][0]

        figtext(0.01, 0.95, r'Minimum Peak$={0:>10.4f}$'.format(best_fit))
        figtext(0.01, 0.93, 'for phi_st, phi_end $({}$, ${})$'.format(min_pa, min_incl))
        figtext(0.01, 0.91, 'at offset $({}$, ${})$'.format(x_offset, y_offset))


        best_fit = np.min(images[1])
        incl_ind, pa_ind = np.where(images[1]==best_fit)
        min_incl, min_pa = yparam[incl_ind][0], xparam[pa_ind][0]

        figtext(0.01, 0.66, r'Min Quadratic Sum$= {0:.4f}$'.format(best_fit))
        figtext(0.01, 0.64, 'for phi_st, phi_end $({}$, ${})$'.format(min_pa, min_incl))

        best_fit = np.min(images[2])
        incl_ind, pa_ind = np.where(images[2]==best_fit)
        min_incl, min_pa = yparam[incl_ind][0], xparam[pa_ind][0]

        figtext(0.01, 0.28, r'Min r>150 AU Quadratic Sum$= {0:.4f}$'.format(best_fit))
        figtext(0.01, 0.26, 'for phi_st, phi_end $({}$, ${})$'.format(min_pa, min_incl))

    plt.tight_layout(h_pad=10.0)

    if save:
        savefig(save)
    if show:
        plt.show()



# Read in image from FITS first moment map
oim = read_mom1_fitsimage(fname='../Plotting/CO_346GHz_line.shifted.robust.5.3sig_mom1.fixed.fits', dpc=dpc, ra0=ra0, dec0=dec0)

# Create square matrix to visualize goodness of fit across models
resid_peak_array = np.zeros((len(yparam),len(xparam)))
peak_offset_array = np.zeros((len(yparam),len(xparam),2))
resid_sum_array = np.zeros((len(yparam),len(xparam)))
resid_sum_out_array = np.zeros((len(yparam),len(xparam)))

# Create model
for cone, aspect, r0, r_st, r_end, v_blob in itertools.product(cones, aspects, r0s, r_sts, r_ends, v_blobs):
    for dpa, incl, grid, pa0, pa_out, i0, i_out, phi_st, phi_end, in itertools.product(dpas, incls, grids, pa0s, pa_outs, i0s, i_outs, phi_sts, phi_ends):

        #Create Model
        images, fit_stats = create_model(image=oim, mstar=mstar, dpc=dpc, incl=incl, dpa=dpa, aspect=aspect, cone=cone, warp=warp, grid=grid, r0=r0, pa0=pa0, pa_out=pa_out, i0=i0, i_out=i_out, blob=blob, phi_st=phi_st, phi_end=phi_end, r_st=r_st, r_end=r_end, v_blob=v_blob)

        # Plot image, model, residuals, and ratio
        model_plot(images=images['image'],
                figsize=[10,15],
                rows=3,
                cols=2,
                extent=[cxmin,cxmax,cymin,cymax],
                x_majmin=[2,0.5],
                y_majmin=[2,0.5],
                fig_titles=['Data', 'Model', 'Residuals', 'Residuals (Velocity Resolution)'],
                xlabel='Relative Right Ascension ()',
                ylabel='Relative Declination ()',
                cbar_each=True,
                cbar_labels=["km/s","km/s","km/s","0.22 km/s increments"],
                interpolation=[B,B,B,N],
                show=True,
                fit_stats=fit_stats,
                save='Shifted_Plots/blob_aspect{}_incl{}_dpa{}_phi_st{}_phi_end{}_r_st{}_r_end{}_v_blob{}.png'.format(aspect, incl, dpa-90., phi_st, phi_end, r_st, r_end, v_blob))
                # save='Shifted_Plots/aspect_{}_cone_{}_dpa_{}_incl{}.png'.format(aspect, cone, dpa-90.0, incl))
                # save='Shifted_Plots/warp_r0_{}_aspect_{}_cone_{}_pa_out_{}_i_out{}_pa0_{}_i0_{}.png'.format(r0, aspect, cone, pa_out-90., i_out, pa0-90., i0))

    # print('FIT')
    # print('FIT')
    # print('FIT')
    # fit_plot(images=[resid_peak_array, resid_sum_array, resid_sum_out_array],
    #         figsize=[10,15],
    #         rows=3,
    #         cols=1,
    #         label_fig=3,
    #         xlabel='Position Angle',
    #         ylabel='Inclination',
    #         x_majmin=[10,1000],
    #         y_majmin=[5,1000],
    #         cpal=cm.Greys,
    #         vmin=[0,0,0],
    #         vmax=[1.3,50,30],
    #         cbar_each=True,
    #         cbar_label=['Peak Residual','Residual Sum', 'Outer Residual Sum'],
    #         fit_stats=peak_offset_array,
    #         show=False,
    #         save='Shifted_Plots/blob_fit_plot_aspect{}_incl{}_dpa{}_r_st{}_r_end{}_v_blob{}.png'.format(aspect,incl, dpa-90., r_st, r_end, v_blob))
            # save='Shifted_Plots/warp_fit_plot_r0_{}_aspect_{}_cone_{}_pa_out_{}_i_out_{}.png'.format(r0, aspect, cone, pa_out-90., i_out))
            # save='Shifted_Plots/fit_plot_aspect_{}_cone_{}.png'.format(aspect, cone))

