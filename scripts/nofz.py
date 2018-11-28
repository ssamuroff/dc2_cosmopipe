import yaml
import numpy as np 
import code.base as base
import sys
import matplotlib.pyplot as plt
import matplotlib
import os
plt.style.use('y1a1')
plt.switch_backend('agg')
matplotlib.rcParams['text.usetex']=False

def pz(z0, sigma=0.05):
	z = np.linspace(0,3.5,1000)
	sig = (1+z0)*sigma
	return np.exp(-(z-z0)*(z-z0)/2/sig/sig)

class generator(base.interface):

	def get_nofz(self, fthin=100):
		# Now let's try a basic mock up for the tomographic binning
		# This comes from Ma, Hu & Huterer (2005) arXiv:0506614
		# But the details aren't really important for the moment

		# Define the redshift binning
		self.edges = self.find_bin_edges(self.info['nbins'])
		ztrue = self.cols['redshift'][self.mask]

		nz = []

		# Cycle through each bin in turn
		for i,(lower,upper) in enumerate(zip(self.edges[:-1],self.edges[1:])):
		    # Select the galaxies with true redshifts in the bin range
		    binsel = (ztrue<upper) & (ztrue>lower)
		    zbin = ztrue[binsel][::fthin]

		    print(i+1,lower,upper,len(zbin))
		    # Now run through the galaxies and stack the per-galaxy PDFs
		    # This is a bad thing to do. I know. Sorry.
		    buffer = np.zeros(1000)

		    for z0 in zbin:
		    	buffer+=pz(z0, sigma=self.info['sigma'])

		    nz.append(buffer)

		self.nz = np.array(nz)

		print('Done')

	def plot(self, show_true=True):
		plt.subplot(1,1,1)
		z = np.linspace(0,3.5,1000)

		for j,n in enumerate(self.nz):
			plt.plot(z,n/np.trapz(n,z),color=self.colours[j])

		for j,(lower,upper) in enumerate(zip(self.edges[:-1], self.edges[1:])):
			plt.axvspan(lower,upper,color=self.colours[j],alpha=0.1)

		if show_true:
			h,b = np.histogram(self.cols['redshift'][self.mask], bins=np.linspace(0,3.5,400), normed=1)
			x = (b[:-1]+b[1:])/2
			plt.plot(x,h,ls=':',color='k')

		plt.xlim(0,2.5)
		plt.ylim(ymin=0)
		plt.yticks(visible=False)
		plt.xticks(visible=True)

		os.system('mkdir -p %s/nofz/'%self.basedir)
		out = np.vstack((z,self.nz))
		np.savetxt('%s/nofz/%s.nz'%(self.basedir,self.sample),out.T)

		plt.subplots_adjust(hspace=0,wspace=0, bottom=0.14)
		plt.xlabel('Redshift $z$', fontsize=16)

		os.system('mkdir -p %s/plots/'%self.basedir)
		plt.savefig('%s/plots/nofz.pdf'%self.basedir)

	def find_bin_edges(self,nbins,w=None):
		"""For an array x, returns the boundaries of nbins equal (possibly weighted by w) bins."""

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
	config = yaml.load(open(sys.argv[-1],'rb'))
	columns = config['basic']['columns'].split()

	nofz = generator(config['basic']['simulation'], columns, basedir=config['basic']['workdir'])
	nofz.create_mask(config)
	nofz.parse_config(config, sections=['nofz'])


	nofz.get_nofz()
	nofz.plot()

main()