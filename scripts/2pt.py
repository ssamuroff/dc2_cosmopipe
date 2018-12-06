import yaml
import numpy as np 
import code.base as base
import sys
import matplotlib.pyplot as plt
import matplotlib
from code.correlation_functions import *
import os
plt.style.use('y1a1')
plt.switch_backend('agg')
matplotlib.rcParams['text.usetex']=False



class correlator(corrfns):
	def __init__(self, cat1, cat2):
		self.corrtype = (cat1.info['ctype'], cat2.info['ctype'])
		print('Correlation type : ', self.corrtype)

		if self.corrtype[0]==self.corrtype[1]:
			self.isauto = True
		else:
			self.isauto = False

		self.bin1, self.p1 = cat1.assign_galaxies_to_bins()
		self.bin2, self.p2 = cat2.assign_galaxies_to_bins()

		self.decide_what_to_do()

		if self.corrtype[0]=='position':
			self.create_random_catalogue(1,cat1)
		if self.corrtype[1]=='position':
			self.create_random_catalogue(2,cat2)

	def process_all(self,cat1,cat2):
		"""Loop through all the possible pairs of tomographic bins
		and compute all of the correlations."""
		done = []
		for i in self.bin1:
			for j in self.bin2:

				# xipm and wtheta are invariant under flipping the 
				# redshift bin indices
				# So we only have to process n(n+1)/2, not n**2 correlations
				if self.isauto and ('%d%d'%(j,i) in done):
					continue

				print('Processing :', i, j)
				self.process(i,j,cat1,cat2)

				done.append('%d%d'%(i,j))

		print('Done.')
		return 0

	def process(self,i,j,cat1,cat2):
		"""Call the chosen function to process a single 
		2pt correlation X(theta), for a single pair of bins i,j """
		if self.corrtype[1]=='position':
			x, creal, cimag, err, err2 = self.workfunction(i,j,cat2,cat1)
		else:
			x, creal, cimag, err, err2 = self.workfunction(i,j,cat1,cat2)
		out = np.array([x,creal,cimag,err,err2])

		os.system('mkdir -p %s/2pt'%(cat1.basedir))
		np.savetxt('%s/2pt/%s_%s_%d_%d.txt'%(cat1.basedir, self.corrtype[0], self.corrtype[1], i, j), out.T)

		return 0


	def create_random_catalogue(self,i,cat):
		"""Generate a catalogue of random points, of the same size
		and drawn from the same area as the object catalogue provided."""

		nrand = cat.cols['ra'].size
		print('Generating random catalogue of %d points'%nrand)

		randoms = np.zeros(nrand, dtype=[('ra',float), ('dec',float)])

		dx = cat.cols['ra'].max() - cat.cols['ra'].min()
		x0 = cat.cols['ra'].min() 
		randoms['ra'] = np.random.rand(nrand) * dx - x0

		dy = cat.cols['dec'].max() - cat.cols['dec'].min()
		y0 = cat.cols['dec'].min() 
		randoms['dec'] = np.random.rand(nrand) * dy - y0

		setattr(self, 'rcat%d'%i, randoms)

		return 0

	def decide_what_to_do(self):
		"""Choose a python function for interfacing with treecorr."""

		s1 = ('shear' in self.corrtype[0]) or ('ellipticity' in self.corrtype[0])
		s2 = ('shear' in self.corrtype[1]) or ('ellipticity' in self.corrtype[1])
		if (s1 and s2):
		   setattr(self, 'workfunction', self.compute_shear_shear)
		   return 0

		p1 = (self.corrtype[0]=='position')
		p2 = (self.corrtype[1]=='position')
		if (p1 and p2):
		   setattr(self, 'workfunction', self.compute_position_position)
		   return 0

		if (s1 and p2) or (p1 and s2):
			setattr(self, 'workfunction', self.compute_position_shear)
			return 0

		print("Unrecognised correlation type. Will assume it's a pair of scalar quantities, matching up to catalogue columns.")
		print("Quite possibly something will break.")
		setattr(self, 'workfunction', self.compute_kappa_kappa)
		return 0

class catalogue(base.interface):

	def assign_galaxies_to_bins(self):
		"""Put galaxies into tomographic bins, according to
		their true redshifts.
		Returns a set of indices ([1,2,3...zbins])
		and an array of length Ngal with the integer bin 
		assignments."""
		self.edges = self.find_bin_edges(self.info['zbins'])
		ztrue = self.cols['redshift']
		zbin = np.zeros(ztrue.size) - 1
		index = []

		# Cycle through each bin in turn
		for i,(lower,upper) in enumerate(zip(self.edges[:-1],self.edges[1:])):
		    # Select the galaxies with true redshifts in the bin range
		    binsel = (ztrue<upper) & (ztrue>lower)
		    zbin[binsel] = i+1
		    index.append(i+1)

		index = np.array(index)

		return index, zbin

	def find_bin_edges(self,nbins,w=None):
		"""For an array x, returns the boundaries of nbins equal 
		(possibly weighted by w) bins."""

		x = self.cols['redshift'][self.mask]

		if w is None:
			xs=np.sort(x)
			r=np.linspace(0.,1.,nbins+1.)*(len(x)-1)
			return xs[r.astype(int)]

		fail=False
		ww=np.sum(w)/nbins
		i=np.argsort(x)
		k=np.linspace(0.,1.,nbins+1.)*(len(x)-1)
		k=k.astype(int)
		r=np.zeros((nbins+1))
		ist=0

		for j in xrange(1,nbins):
			if k[j]<r[j-1]:
				print('Random weight approx. failed - attempting brute force approach')
				fail=True
				break

			w0=np.sum(w[i[ist:k[j]]])
			if w0<=ww:
				for l in xrange(k[j],len(x)):
					w0+=w[i[l]]
					if w0>ww:
						r[j]=x[i[l]]
						ist=l
						break
			else:
				for l in xrange(k[j],0,-1):
					w0-=w[i[l]]
					if w0<ww:
						r[j]=x[i[l]]
						ist=l
						break

		if fail:
			ist=np.zeros((nbins+1))
			ist[0]=0
			for j in xrange(1,nbins):
				wsum=0.
				for k in xrange(ist[j-1].astype(int),len(x)):
					wsum+=w[i[k]]
					if wsum>ww:
						r[j]=x[i[k-1]]
						ist[j]=k
						break

		r[0]=x[i[0]]
		r[-1]=x[i[-1]]
		return r




def main():
	# Two config files defining two samples to correlate
	config1 = yaml.load(open(sys.argv[-2],'rb'))
	config2 = yaml.load(open(sys.argv[-1],'rb'))

	columns1 = config1['basic']['columns'].split()
	columns2 = config2['basic']['columns'].split()

	cat1 = catalogue(config1['basic']['simulation'], columns1, basedir=config1['basic']['workdir'])
	cat1.parse_config(config1, sections=['nofz','2pt'])
	cat1.create_mask(config1)

	cat2 = catalogue(config2['basic']['simulation'], columns2, basedir=config2['basic']['workdir'])
	cat2.parse_config(config2, sections=['nofz','2pt'])
	cat2.create_mask(config2)

	# Set up an object to handle the correlations
	corr = correlator(cat1, cat2)

	corr.process_all(cat1,cat2)


main()