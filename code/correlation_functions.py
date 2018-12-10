import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
plt.style.use('y1a1')
matplotlib.rcParams['text.usetex']=False
import treecorr

colnames = {'shear' : ('shear_1', 'shear_2'), 'ellipticity' : ('ellipticity_1', 'ellipticity_2'), 'ellipticity_true' : ('ellipticity_1_true', 'ellipticity_2_true')}

class corrfns:
	def compute_shear_shear(self, i, j, cat1, cat2):
		maski = (cat1.mask) & (self.p1==i)
		maskj = (cat2.mask) & (self.p2==j)

		# Initialised the catlogues
		namei_1,namei_2 = colnames[self.corrtype[0]]
		cat_i = treecorr.Catalog(g1=cat1.cols[namei_1][maski], g2=cat1.cols[namei_2][maski], ra=cat1.cols['ra'][maski], dec=cat1.cols['dec'][maski], ra_units='deg', dec_units='deg')

		namej_1,namej_2 = colnames[self.corrtype[1]]
		cat_j = treecorr.Catalog(g1=cat2.cols[namej_1][maskj], g2=cat2.cols[namej_2][maskj], ra=cat2.cols['ra'][maskj], dec=cat2.cols['dec'][maskj], ra_units='deg', dec_units='deg')

		# Set up the correlation
		# Note that we're using the binning configuration from 
		# the first of the two config files
		# might be good to check and give a warning if the two files
		# specify different 2pt binning parameters
		gg = treecorr.GGCorrelation(nbins=cat1.info['tbins'], min_sep=cat1.info['tmin'], max_sep=cat1.info['tmax'], sep_units='arcmin', bin_slop=0.1, verbose=True,num_threads=1)

        # And process it
		gg.process(cat_i,cat_j)
        
		theta = np.exp(gg.meanlogr)
		xip = gg.xip
		xim = gg.xim
		xiperr = ximerr = np.sqrt(gg.varxi)

		return theta, xip, xim, xiperr, ximerr

	def compute_position_shear(self, i, j, cat1, cat2):
		maski = (cat1.mask) & (self.p1==i)
		maskj = (cat2.mask) & (self.p2==j)

		# Initialised the catlogues
		cat_i = treecorr.Catalog(ra=cat1.cols['ra'][maski], dec=cat1.cols['dec'][maski], ra_units='deg', dec_units='deg')

		rcat_i = treecorr.Catalog(ra=self.rcat1['ra'][maski], dec=self.rcat1['dec'][maski], ra_units='deg', dec_units='deg')

		namej_1,namej_2 = colnames[self.corrtype[1]]
		cat_j = treecorr.Catalog(g1=cat2.cols[namej_1][maskj], g2=cat2.cols[namej_2][maskj], ra=cat2.cols['ra'][maskj], dec=cat2.cols['dec'][maskj], ra_units='deg', dec_units='deg')

		# Set up the correlation
		ng = treecorr.NGCorrelation(nbins=cat1.info['tbins'], min_sep=cat1.info['tmin'], max_sep=cat1.info['tmax'], sep_units='arcmin', bin_slop=0.1, verbose=True,num_threads=1)
		rg = treecorr.NGCorrelation(nbins=cat1.info['tbins'], min_sep=cat1.info['tmin'], max_sep=cat1.info['tmax'], sep_units='arcmin', bin_slop=0.1, verbose=True,num_threads=1)

        # And process it
		ng.process(cat_i,cat_j)
		rg.process(rcat_i,cat_j)

		gammat, gammat_im, gammaterr = ng.calculateXi(rg)

		theta = np.exp(ng.meanlogr)
		gammaterr = np.sqrt(gammaterr)

		return theta, gammat, gammat_im, gammaterr, gammaterr

	def compute_position_position(self, i, j, cat1, cat2):
		maski = (cat1.mask) & (self.p1==i)
		maskj = (cat2.mask) & (self.p2==j)
		rmaski = np.random.choice(self.rcat1['ra'].size, size=maski[maski].size, replace=False)
		rmaskj = np.random.choice(self.rcat2['ra'].size, size=maskj[maskj].size, replace=False)

		# Initialised the catlogues
		cat_i = treecorr.Catalog(ra=cat1.cols['ra'][maski], dec=cat1.cols['dec'][maski], ra_units='deg', dec_units='deg')

		cat_j = treecorr.Catalog(ra=cat2.cols['ra'][maskj], dec=cat2.cols['dec'][maskj], ra_units='deg', dec_units='deg')

		rancat_i = treecorr.Catalog(ra=self.rcat1['ra'][rmaski], dec=self.rcat1['dec'][rmaski], ra_units='deg', dec_units='deg')
		rancat_j = treecorr.Catalog(ra=self.rcat2['ra'][rmaskj], dec=self.rcat2['dec'][rmaskj], ra_units='deg', dec_units='deg')

		# Trigger a warning if the random catalogues differ significantly from
		# the main catalogues in size 
		checki = (abs(cat_i.x.size - rancat_i.x.size) * 1./cat_i.x.size)>0.25
		checkj = (abs(cat_j.x.size - rancat_j.x.size) * 1./cat_j.x.size)>0.25

		if checki or checkj:
			print("Warning: there are either significantly more or fewer randoms than actual galaxies in one or both samples.")

		# Set up the correlation
		nn = treecorr.NNCorrelation(nbins=cat1.info['tbins'], min_sep=cat1.info['tmin'], max_sep=cat1.info['tmax'], sep_units='arcmin', bin_slop=0.1, verbose=True,num_threads=1)

		# And process it
		nn.process(cat_i,cat_j)

		nr = treecorr.NNCorrelation(nbins=cat1.info['tbins'], min_sep=cat1.info['tmin'], max_sep=cat1.info['tmax'], sep_units='arcmin', bin_slop=0.1, verbose=True,num_threads=1)
		rn = treecorr.NNCorrelation(nbins=cat1.info['tbins'], min_sep=cat1.info['tmin'], max_sep=cat1.info['tmax'], sep_units='arcmin', bin_slop=0.1, verbose=True,num_threads=1)
		rr = treecorr.NNCorrelation(nbins=cat1.info['tbins'], min_sep=cat1.info['tmin'], max_sep=cat1.info['tmax'], sep_units='arcmin', bin_slop=0.1, verbose=True,num_threads=1)

		nr.process(cat_i,rancat_j)
		rn.process(rancat_i,cat_j)
		rr.process(rancat_i,rancat_j)

		wtheta,wthetaerr = nn.calculateXi(rr,dr=nr,rd=rn)
		theta = np.exp(nn.meanlogr)
		wthetaerr = np.sqrt(wthetaerr)

		return theta, wtheta, np.array([0]*len(theta)), wthetaerr, np.array([0]*len(theta))
