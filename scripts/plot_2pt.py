import yaml
import numpy as np 
import code.base as base
import sys
import matplotlib.pyplot as plt
import matplotlib
import os
import glob
plt.style.use('y1a1')
plt.switch_backend('pdf')
matplotlib.rcParams['text.usetex']=False
matplotlib.rcParams["ytick.minor.visible"]=False
matplotlib.rcParams["ytick.minor.width"]=0.1

class cornerplot:
    def __init__(self, options):
        self.basedir = options['basic']['workdir']
        files = glob.glob('%s/2pt/shear_shear_*txt'%self.basedir)

        # Work out how many bins we have from the file names  
        # Sorry, this is a bit fiddly...
        self.n1 = np.unique([int(os.path.basename(f).split('_')[2]) for f in files]).size
        self.n2 = np.unique([int(os.path.basename(f).split('_')[3].replace('.txt','')) for f in files]).size

    def get_spectra(self, i, j):
        data = np.loadtxt('%s/2pt/shear_shear_%d_%d.txt'%(self.basedir, i+1, j+1)).T
        x = data[0] # should already be in arcminutes
        xip = data[1]
        dxip = data[3]
        xim = data[2]
        dxim = data[4]
        return x, xip, xim, dxip, dxim

    def make(self):

        rows, cols = self.n1, self.n2
        count = 0

        for i in range(self.n1):
            for j in range(self.n2):
                count+=1
                if j>i:
                    continue

                print(i,j)
                x, xip,xim, dp, dm = self.get_spectra(i,j)

                #posp = positions[(i+1,j+1,"+")]
                ax = plt.subplot(rows,cols,count)

                ax.annotate("(%d, %d)"%(i+1,j+1), (3,4.9), textcoords='data', fontsize=9 )
                ax.yaxis.set_tick_params(which='minor', left='off', right='off')
                plt.ylim(-0.5,6.8)
                #plt.yscale("log")
                plt.xscale("log")
                plt.yticks(fontsize=9)

                if j==self.n2-1:
                    #plt.ylabel(r"$\theta \xi_+(\theta)$ / arcmin", fontsize=12)
                    plt.xlabel(r"$\theta$ / arcmin", fontsize=8)
                    plt.yticks(fontsize=9)

                if j==0:
                    plt.yticks([0,2,4,6],['0', '2', '4', '6'])
                    if i==self.n1-2:
                        plt.ylabel(r"$\theta \xi_\pm(\theta) / \times 10^{-4}$ arcmin", fontsize=9)

                else:
                    plt.yticks(visible=False)

                plt.xlim(2.2,270)
                plt.xticks([10,100],["10", "100"], fontsize=9)
                
                plt.axhline(0, color='k', ls=':')

                plt.errorbar(x, 1e4 * x * xip, yerr=dp, ls="none", marker="D", markersize=4,  ecolor="purple", markeredgecolor="purple", markerfacecolor="none", label=r'$\xi_+$')
                plt.errorbar(x, 1e4 * x * xim, yerr=dm, ls="none", marker="o", markersize=4, ecolor="pink", markeredgecolor="pink", markerfacecolor="pink", label=r'$\xi_-$')

                if (i==0) and (j==0):
                    plt.legend(bbox_to_anchor=(1.8, 1), fontsize=9)


        plt.subplots_adjust(hspace=0,wspace=0.0, bottom=0.14, top=0.97, left=0.14)
        plt.savefig("%s/plots/cornerplot-shear_shear.pdf"%self.basedir)

def main():
    config = yaml.load(open(sys.argv[-1],'rb'))
    plotter = cornerplot(config)
    plotter.make()

main()