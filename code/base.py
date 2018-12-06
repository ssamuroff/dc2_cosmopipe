import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
plt.style.use('y1a1')
matplotlib.rcParams['text.usetex']=False

import GCRCatalogs
from GCR import GCRQuery



# Start by defining a class to handle interactions with the external DC2 code/data
class interface:
	def __init__(self, simulation, columns, basedir=0):
		self.colours=['purple', 'pink', 'plum', 'royalblue', 'forestgreen', 'steelblue', 'orange', 'k']
		self.simulation = simulation

		if isinstance(basedir,int):
			self.basedir = '/global/cscratch1/sd/sws/dc2/data/%s/desy5/'%simulation
		else:
			self.basedir = basedir

		print("Using %s as working directory."%self.basedir)
		print("If it didn't exist before, it does now.")
		os.system('mkdir -p %s'%self.basedir)

		print('Loading catalogue')
		cat = GCRCatalogs.load_catalog(simulation)

		print('Retrieving (%d) requested column(s).'%len(columns))
		self.cols = cat.get_quantities(columns)

		print('Catalogue contains %d M objects'%(len(self.cols[ columns[0] ])*1.0/1e6 ))
		self.ngal = len(self.cols[ columns[0] ])

	def create_mask(self, options):
		"""Generates a selection mask to apply to the catalogues.
		   The cuts must be provided in the config file in the form 
		   column_name1,value1,+/- column_name2,value2,+/-
		   (the spacing and commas are important) """

		cuts = options['basic']['cuts'].split()

		self.mask = np.ones(self.ngal).astype(bool)

		for cut in cuts:
			colname,value,direction = cut.split(',')

			value = float(value)

			if (direction=='+'):
				mask = (self.cols[colname]>value) & np.isfinite(self.cols[colname])
			elif (direction=='-'):
				mask = (self.cols[colname]<value) & np.isfinite(self.cols[colname])

			self.mask = self.mask & mask

		if ('colour_split' in options['basic'].keys()) and (not options['basic']['colour_split'] is None):
			csplit = options['basic']['colour_split']
			
			self.true_colour_split(csplit)

		print( 'Mask leaves %d/%d galaxies'%( len(self.cols['redshift'][self.mask]), len(self.cols['redshift'])) )

		return 0

	def colour_split(self, colour, show=True):
		print('Colour bin : %s'%colour)

		if show:
			fig = plt.figure(figsize=(9, 9), dpi= 80, facecolor='w', edgecolor='k')

		# Until I can come up with anything better, I'm just going to hard code these numbers
		# to the values derived here https://arxiv.org/abs/1811.06989
		a= [0.037,0.12,0.05,0.0]
		c0=[-0.1,-1.7,0.15,1.6]

		redflag = np.zeros(self.cols['redshift'][self.mask].size)

		rmag = self.cols['mag_r_lsst'][self.mask]
		zmag = self.cols['mag_z_lsst'][self.mask]
		rz = rmag - zmag

		ybins = np.linspace(0.3,2.5,70)
		xbins = np.linspace(20,25.5,70)

		self.edges = self.find_bin_edges(self.info['zbins'])
		ztrue = self.cols['redshift'][self.mask]

		# Cycle through each bin in turn
		for i,(lower,upper) in enumerate(zip(self.edges[:-1],self.edges[1:])):

			if i+1!=self.nbins:
				plt.xticks(visible=False)

			mask = (ztrue<upper) & (ztrue>lower)

			x0 = rmag[mask]
			y0 = rz[mask]

			if show:
				plt.subplot('%d1%d'%(self.nbins,i+1), aspect=0.85)

			if show:
				xlim = (18,20.0)
				ylim = (0.3,1.2)

				counts,xbins,ybins, = np.histogram2d(x0, y0, bins=50, range=[xlim,ylim], normed=1 )
				y=(ybins[:-1]+ybins[1:])/2
				x=(xbins[:-1]+xbins[1:])/2
				xx,yy=np.meshgrid(x,y)

				print("Making histograms...")
				C = plt.contour(xx,yy,counts.T, 4, colors='purple',linestyles='-', linewidth=.5)
				plt.xlim(xbins.min(),xbins.max())
				plt.ylim(ybins.min(),ybins.max())
				#plt.yticks([0.5,1.0,1.5,2.0])

				xl = np.linspace(xlim[0],xlim[1],100)
				#lin = a[i]*xl+c0[i]

				#plt.plot(xl,lin, color="forestgreen", ls="--", lw=1.5)

				plt.annotate('(%d)'%(i+1), xy=(20.2,2), fontsize=18)
				plt.ylabel('$r-z$')

			#cmask = y0>(x0*a[i]+c0[i])
			#blank = np.zeros(redflag[mask].size)
			#blank[cmask]=1
			#redflag[mask] = blank

			#print(i+1, '%3.3f M/ %3.3f M = %3.3f'%(len(blank[cmask])/1e6,redflag[mask].size/1e6,len(blank[cmask])*1/redflag[mask].size))


		if show:
			plt.subplots_adjust(wspace=0, hspace=0)
			plt.xlabel('$r$-band Magnitude')
			plt.savefig('%s/plots/colour_mag.pdf'%self.basedir)

		fred = len(redflag[redflag==1])*1./len(redflag)
		print('')
		print('Global red fraction : %3.3f'%fred)

		newmask = np.zeros(self.mask.size).astype(bool)
		newmask[self.mask] = redflag

		print('Updating mask.')
		if colour=='red':
			self.mask = self.mask & newmask
		elif colour=='blue':
			self.mask = self.mask & np.invert(newmask)

		return 0



	def parse_config(self, config, sections=None):

		self.sample = config['basic']['sample']
		self.nbins = config['nofz']['zbins']

		self.info = {}

		if sections is None:
			sections = config.keys()

		for section in sections:
			for name in config[section].keys():
				self.info[name] = config[section][name]
				print('%s = '%name, self.info[name])

		return None

	def true_colour_split(self, colour, show=True):
		print('Colour bin : %s'%colour)

		if show:
			fig = plt.figure(figsize=(9, 9), dpi= 80, facecolor='w', edgecolor='k')


		#redflag = np.zeros(self.cols['redshift'][self.mask].size)

		rmag = self.cols['mag_true_r_lsst'][self.mask]
		zmag = self.cols['mag_true_z_lsst'][self.mask]
		rz = rmag - zmag

		ybins = np.linspace(-0.5,2.5,70)
		xbins = np.linspace(20,25.5,70)

		a0 = 0.037
		c0 = -0.175
		lineatx = (self.cols['mag_true_r_lsst']*a0 + c0)
		rz0 = self.cols['mag_true_r_lsst'] - self.cols['mag_true_z_lsst']

		redflag = rz0 < lineatx
		fred = rz0[self.mask & redflag].size * 1.0 / rz0[self.mask].size

		print('Global red fraction: %3.3f percent (%d/%d)'%(fred * 100., rz0[self.mask & redflag].size , rz0[self.mask].size ))

		if colour=='red':
			self.mask = self.mask & redflag
		elif colour=='blue':
			self.mask = self.mask & np.invert(redflag)

		if show:
				xlim = (18,25.5)
				ylim = (-0.5,2.5)

				counts,xbins,ybins = np.histogram2d(rmag, rz, bins=50, range=[xlim,ylim], normed=1 )
				h1d, b = np.histogram(rz, bins=np.linspace(ylim[0],ylim[1],60), normed=1 )
				y=(ybins[:-1]+ybins[1:])/2
				x=(xbins[:-1]+xbins[1:])/2
				xx,yy=np.meshgrid(x,y)

				print("Making histograms...")

				plt.close()
				plt.plot((b[1:]+b[:-1])/2,h1d, color='purple', ls='-', lw=2)
				plt.xlabel('$r-z$', fontsize=16)
				plt.yticks(visible=False)
				plt.xlim(ylim[0],ylim[1])
				plt.ylim(ymin=0)
				plt.subplots_adjust(bottom=0.15)
				plt.savefig('%s/plots/rz_colour_hist_true.pdf'%self.basedir)


				plt.close()
				C = plt.contour(xx,yy,counts.T, 8, colors='purple',linestyles='-', linewidth=.5)
				plt.xlim(xbins.min(),xbins.max())
				plt.ylim(ybins.min(),ybins.max())
				plt.ylabel('$r-z$', fontsize=16)
				plt.xlabel('$r$', fontsize=16)
				plt.subplots_adjust(bottom=0.15, left=0.15)

				xl = np.linspace(xlim[0],xlim[1],100)
				lin = a0*xl+c0

				plt.plot(xl,lin, color="forestgreen", ls="--", lw=1.5)
				plt.savefig('%s/plots/rz_rmag_2dhist_true.pdf'%self.basedir)
				#plt.yticks([0.5,1.0,1.5,2.0])
				plt.close()

		return 0



def show_tables():
	"""Print out a numbered list of catalogues available."""

	tables = GCRCatalogs.available_catalogs

	for i,name in tables.keys():
		print(i,name)
	return 0