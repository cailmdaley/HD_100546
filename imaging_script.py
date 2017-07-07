# Cail Daley
# My iteration of Catherine Walsh's imaging script, for use on HD 100546.
#=============================================================================#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Comparison of different start velocity clean beams

# Original Velocity, Original Briggs CLEAN (robust=0.5):
# restoring beam: 0.93555 X 0.42365 arcsec

# Original Velocity, Natural CLEAN:
# restoring beam: 0.99896 X 0.50444 arcsec

# Original Velocity, Uniform CLEAN:
# restoring beam: 0.9236 X 0.3811 arcsec



# ALL of these have robust=0
# Negative Shifted Velocity, Briggs CLEAN (-0.075):
# restoring beam: 0.92176 X 0.39266 arcsec

# Original Velocity, Briggs CLEAN, redone:
# restoring beam: 0.92176 X 0.39266 arcsec

# Shifted Velocity, Briggs CLEAN (+0.075):
# restoring beam: 0.92175 X 0.39264 arcsec

# Shifted More Velocity, Briggs CLEAN (+0.15):
# restoring beam: 0.92176 X 0.39266 arcsec
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# make a dirty image to estimate rms noise in emission free channels?
# upon visual inspection there is significant emission in +12 km/s channel
# source velocity is ~ 5.5 km/s
# shift velocity axis to -6 km/s --> +18 km/s
clean(vis='sb2_spw3_corrected.fixed.ms.contsub',
	imagename='CO_346GHz_line.dirty',
	imagermode='csclean',
	mode='velocity',
	start='-5.925km/s',
	width='0.15km/s',
	nchan=160,
	outframe='LSRK',
	restfreq='345.7959899GHz',
	interactive=F,
	imsize=[256,256],
	cell='0.12arcsec',
	weighting='natural',
	niter=0)


#box width X height: 7.18012 X 7.63456
#box center x X y:  11:33:25.310 X -70.11.41.788
#rms estimated within box in line free channels of dirty map from -6 km/s to -2.4 km/s (channel 0-25); source not quite detectable at this velocity but rms steadily climbs from this point on.
#25 channels:
#1.542479e-02+1.539211e-02+1.750176e-02+1.603433e-02+1.679388e-02+1.679546e-02+1.786114e-02+1.827202e-02+1.577174e-02+1.728114e-02+1.474196e-02+1.665305e-02+1.703183e-02+1.832825e-02+1.742577e-02+1.550170e-02+1.633467e-02+1.591223e-02+1.697805e-02+1.720553e-02+1.705352e-02+1.968431e-02+1.684617e-02+1.787351e-02+1.857922e-02

#rms estimated within box in line free channels from 17.85 km/s ---> 14.55 km/s (channel 137-160); source not quite visually detectable at this velocity but rms steadily climbs from this point on.
#23 channels:
#1.522460e-02+1.854534e-02+1.400986e-02+1.610523e-02+1.713309e-02+1.605422e-02+1.757498e-02+1.733189e-02+1.627310e-02+1.596563e-02+1.737688e-02+1.728297e-02+1.722569e-02+1.631972e-02+1.669923e-02+1.783142e-02+1.859657e-02+1.733985e-02+1.792821e-02+1.725908e-02+1.777731e-02+1.825241e-02+1.698060e-02

#Average rms = (0.42327813999999997 + 0.39108788) / (25 + 23) = 0.01696595875 Jy = 16.97 mJy ~ 17 mJy



#BRIGGS CLEAN
#interactive clean with channel-independent mask:
clean(vis='sb2_spw3_corrected.fixed.ms.contsub',
	imagename='CO_346GHz_line.independent1',
	imagermode='csclean',
	mode='velocity',
	start='-6km/s',
	width='0.15km/s',
	nchan=220,
	outframe='LSRK',
	restfreq='345.7959899GHz',
	interactive=T,
    mask='CO_346GHz_line_full.mask',
	imsize=[256,256],
	cell='0.12arcsec',
	weighting='briggs',
	robust=0.5,
	niter=20000,
    threshold = '1mJy')
imstat(imagename = 'CO_346GHz_line.independent1.image/',
    box = '0,0,255,90,166,0,255,255',
    chans = '25~136')
# max residual = .062 Jy
# peak flux density = 3.06961 Jy/beam
# total flux = 966.621 Jy
# rms = 18.5 mJy ~ 19 mJy (about the same as Catherine's notes)
# S/N = 166 (compare with Catherine's value of 162)
# restoring beam: 0.9356 X 0.4237 arcsec


# Make 0th moment with no clip, find rms outside source to create a 3 sigma clip
immoments(imagename='CO_346GHz_line.independent1.image',
    moments=[0,1,2],
    outfile='CO_346GHz_line.independent1.uncut_mom',
    chans='25~136')
#imstat with same box as before to get rms:
imstat(imagename = 'CO_346GHz_line.independent1.uncut_mom0/',
    box = '0,0,255,90,166,0,255,255')
# rms = 0.054 Jy/beam * km/s
# Make region defined by 3 sigma contour of the uncut 0th moment map:
# Region = CO_346GHz_line.independent1.3sig_region

# Make disk-integrated line profile, export to fits file:
#Profile = CO_346GHz_line.independent1.spectral_profile.fits


#create moment maps with 3 sigma clip
immoments(imagename='CO_346GHz_line.independent1.image',
    moments=[0,1,2],
    outfile='CO_346GHz_line.independent1.3sig_mom',
    chans='25~136',
    includepix = [3 * 0.054,1000])

#scale velocities so that the source velocity lies at 0 km/s:
immath(imagename='CO_346GHz_line.independent1.3sig_mom1',
        mode = 'evalexpr',
        expr = 'IM0 - 5.7',
        outfile = 'CO_346GHz_line.independent1.3sig_mom1.fixed')

# Create major axis position-velocity slice:
# Position angle of HD 100546 is 146 degrees
# Image size is 256 x 256 pixels

# Slice along disk major axis needs to be at least 2 pixels away from image edges
# x = [64,192], y = [33,223]
mypv = impv(imagename ='CO_346GHz_line.independent1.image',
    outfile = 'CO_346GHz_line.independent1.pv_major.im',
	start=[64,33], end=[192,223],
	overwrite=T)






#UNIFORM CLEAN
#interactive clean with channel-independent mask:
clean(vis='sb2_spw3_corrected.fixed.ms.contsub',
	imagename='CO_346GHz_line.independent_uniform',
	imagermode='csclean',
	mode='velocity',
	start='-6km/s',
	width='0.15km/s',
	nchan=220,
	outframe='LSRK',
	restfreq='345.7959899GHz',
	interactive=T,
    mask='CO_346GHz_line_full.mask',
	imsize=[256,256],
	cell='0.12arcsec',
	weighting='uniform',
	niter=20000,
    threshold = '1mJy')
imstat(imagename = 'CO_346GHz_line.independent_uniform.image',
    box = '0,0,255,90,166,0,255,255',
    chans = '25~136')
# max residual = .066 Jy
# peak flux density = 2.78734 Jy/beam
# total flux = 926.891 Jy
# rms = 27.71 mJy ~ 28 mJy (3/2 as large as with briggs clean)
# S/N = 101 (compare with Catherine's value of 162)
# restoring beam: 0.9236 X 0.3811 arcsec (compare with briggs beam of 0.9356 X 0.4237 arcsec)


# Make 0th moment with no clip, find rms outside source to create a 3 sigma clip
immoments(imagename = 'CO_346GHz_line.independent_uniform.image',
    moments=[0,1,2],
    outfile='CO_346GHz_line.independent_uniform.uncut_mom',
    chans='25~136')
#imstat with same box as before to get rms:
imstat(imagename = 'CO_346GHz_line.independent_uniform.uncut_mom0/',
    box = '0,0,255,90,166,0,255,255')
# rms = 0.072 Jy/beam * km/s
# Make region defined by 3 sigma contour of the uncut 0th moment map:
# Region = CO_346GHz_line.independent_uniform.3sig_region

# Make disk-integrated line profile, export to fits file:
#Profile = CO_346GHz_line.independent_uniform.spectral_profile.fits


#create moment maps with 3 sigma clip
immoments(imagename='CO_346GHz_line.independent_uniform.image',
    moments=[0,1,2],
    outfile='CO_346GHz_line.independent_uniform.3sig_mom',
    chans='25~136',
    includepix = [3 * 0.072,1000])

#scale velocities so that the source velocity lies at 0 km/s:
immath(imagename='CO_346GHz_line.independent_uniform.3sig_mom1',
        mode = 'evalexpr',
        expr = 'IM0 - 5.7',
        outfile = 'CO_346GHz_line.independent_uniform.3sig_mom1.fixed')


# Create major axis position-velocity slice:
# Position angle of HD 100546 is 146 degrees
# Image size is 256 x 256 pixels

# Slice along disk major axis needs to be at least 2 pixels away from image edges
# x = [64,192], y = [33,223]
mypv = impv(imagename ='CO_346GHz_line.independent_uniform.image',
    outfile = 'CO_346GHz_line.independent_uniform.pv_major.im',
	start=[64,33], end=[192,223],
	overwrite=T)



#Briggs CLEAN robust=0, shifted +0.075 start velocity
#interactive clean with channel-independent mask:
clean(vis='sb2_spw3_corrected.fixed.ms.contsub',
	imagename='CO_346GHz_line.shifted',
	imagermode='csclean',
	mode='velocity',
	start='-5.925km/s',
	width='0.15km/s',
	nchan=220,
	outframe='LSRK',
	restfreq='345.7959899GHz',
	interactive=T,
    mask='CO_346GHz_line_full.mask',
	imsize=[256,256],
	cell='0.12arcsec',
	weighting='briggs',
	niter=20000,
    threshold = '1mJy')
imstat(imagename = 'CO_346GHz_line.independent_shifted.image/',
    box = '0,0,255,90,166,0,255,255',
    chans = '25~136')
# max residual = .061 Jy
# peak flux density = 2.855 Jy/beam
# total flux = 942.752 Jy
# rms = 21.8 mJy ~ 22 mJy (a few mJy higher than Catherine's notes)
# S/N = 131 (compare with Catherine's value of 162)
# restoring beam: 0.92175 X 0.39264 arcsec


# Make 0th moment with no clip, find rms outside source to create a 3 sigma clip
immoments(imagename='CO_346GHz_line.independent_shifted.image',
    moments=[0,1,2],
    outfile='CO_346GHz_line.independent_shifted.uncut_mom',
    chans='25~136')
#imstat with same box as before to get rms:
imstat(imagename = 'CO_346GHz_line.independent1.uncut_mom0/',
    box = '0,0,255,90,166,0,255,255')
# rms = 0.060 Jy/beam * km/s
# Make region defined by 3 sigma contour of the uncut 0th moment map:
# Region = CO_346GHz_line.independent_shifted.3sig_region

# Make disk-integrated line profile, export to fits file:
#Profile = CO_346GHz_line.independent1.spectral_profile.fits


#create moment maps with 3 sigma clip
immoments(imagename='CO_346GHz_line.independent_shifted.image',
    moments=[0,1,2],
    outfile='CO_346GHz_line.independent_shifted.3sig_mom',
    chans='25~136',
    includepix = [3 * 0.06,1000])

#scale velocities so that the source velocity lies at 0 km/s:
immath(imagename='CO_346GHz_line.independent_shifted.3sig_mom1',
        mode = 'evalexpr',
        expr = 'IM0 - 5.625',
        outfile = 'CO_346GHz_line.independent_shifted.3sig_mom1.fixed')




#Briggs CLEAN robust=0.5, shifted +0.075 start velocity
clean(vis='sb2_spw3_corrected.fixed.ms.contsub',
	imagename='CO_346GHz_line.shifted.robust.5',
	imagermode='csclean',
	mode='velocity',
	start='-5.925km/s',
	width='0.15km/s',
	nchan=220,
	outframe='LSRK',
	restfreq='345.7959899GHz',
	interactive=T,
    mask='CO_346GHz_line_full.mask',
	imsize=[256,256],
	cell='0.12arcsec',
	weighting='briggs',
	robust=0.5,
	niter=20000,
    threshold = '1mJy')
imstat(imagename = 'CO_346GHz_line.shifted.robust.5.image/',
    box = '0,0,255,90,166,0,255,255',
    chans = '25~136')
# max residual = .060 Jy
# peak flux density = 3.06271 Jy/beam
# total flux = 963.563 Jy
# rms = 18.5 mJy ~ 19 mJy (about the same as Catherine's notes)
# S/N = 166 (compare with Catherine's value of 162)
# restoring beam: 0.9356 X 0.4237 arcsec


# Make 0th moment with no clip, find rms outside source to create a 3 sigma clip
immoments(imagename='CO_346GHz_line.shifted.robust.5.image',
    moments=[0,1,2],
    outfile='CO_346GHz_line.shifted.robust.5.uncut_mom',
    chans='25~136')
#imstat with same box as before to get rms:
imstat(imagename = 'CO_346GHz_line.shifted.robust.5.uncut_mom0/',
    box = '0,0,255,90,166,0,255,255')
# rms = 0.055 Jy/beam * km/s
# Make region defined by 3 sigma contour of the uncut 0th moment map:
# Region = CO_346GHz_line.shifted.robust.5.3sig_region

# Make disk-integrated line profile, export to fits file:
#Profile = CO_346GHz_line.shifted.robust.5.spectral_profile.fits


#create moment maps with 3 sigma clip
immoments(imagename='CO_346GHz_line.shifted.robust.5.image',
    moments=[0,1,2],
    outfile='CO_346GHz_line.shifted.robust.5.3sig_mom',
    chans='25~136',
    includepix = [3 * 0.055,1000])

#scale velocities so that the source velocity lies at 0 km/s:
immath(imagename='CO_346GHz_line.shifted.robust.5.3sig_mom1',
        mode = 'evalexpr',
        expr = 'IM0 - 5.625',
        outfile = 'CO_346GHz_line.shifted.robust.5.3sig_mom1.fixed')

# Create major axis position-velocity slice:
# Position angle of HD 100546 is 146 degrees
# Image size is 256 x 256 pixels

# Slice along disk major axis needs to be at least 2 pixels away from image edges
# x = [59,197], y = [33,223]
mypv = impv(imagename ='CO_346GHz_line.shifted.robust.5.image',
    outfile = 'CO_346GHz_line.shifted.robust.5.pv_major.im',
	start=[59,33], end=[197,223],
	overwrite=T)
