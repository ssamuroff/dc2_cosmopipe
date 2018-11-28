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

		print('Retrieving (%d) requested columns.'%len(columns))
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
				mask = (self.cols[colname]>value)
			elif (direction=='-'):
				mask = (self.cols[colname]<value)

			self.mask = self.mask & mask

		print( 'Mask leaves %d/%d galaxies'%( len(self.cols['redshift'][self.mask]), len(self.cols['redshift'])) )

		return 0


	def parse_config(self, config, sections=None):

		self.sample = config['basic']['sample']

		self.info = {}

		if sections is None:
			sections = config.keys()

		for section in sections:
			for name in config[section].keys():
				self.info[name] = config[section][name]
				print('%s = '%name, self.info[name])

		return None



def show_tables():
	"""Print out a numbered list of catalogues available."""

	tables = GCRCatalogs.available_catalogs

	for i,name in tables.keys():
		print(i,name)
	return 0