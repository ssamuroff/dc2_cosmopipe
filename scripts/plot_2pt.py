import yaml
import numpy as np 
import code.base as base
import sys
import matplotlib.pyplot as plt
import matplotlib
import os
import glob
import argparse
plt.style.use('y1a1')
plt.switch_backend('pdf')
matplotlib.rcParams['text.usetex']=False
matplotlib.rcParams["ytick.minor.visible"]=False
matplotlib.rcParams["ytick.minor.width"]=0.1

class cornerplot:
    def __init__(self, options, args):
        self.basedir = options['basic']['workdir']
        self.files_xipm = glob.glob('%s/2pt/shear_shear_*txt'%self.basedir)
        self.files_gammat = glob.glob('%s/2pt/position_shear_*txt'%self.basedir)
        self.files_wtheta = glob.glob('%s/2pt/position_position_*txt'%self.basedir)

        # Work out how many bins we have from the file names  
        # Sorry, this is a bit fiddly...
        self.nsrc = np.unique([int(os.path.basename(f).split('_')[2]) for f in self.files_xipm]).size
        self.nlens = np.unique([int(os.path.basename(f).split('_')[2]) for f in self.files_wtheta]).size
        #self.n2 = np.unique([int(os.path.basename(f).split('_')[3].replace('.txt','')) for f in files]).size

        self.ss = args.xipm
        self.ns = args.gammat
        self.nn = args.wtheta

    def get_spectra(self, i, j, ptype='shear_shear'):
        data = np.loadtxt('%s/2pt/%s_%d_%d.txt'%(self.basedir, ptype, i+1, j+1)).T
        x = data[0] # should already be in arcminutes
        xip = data[1]
        dxip = data[3]
        xim = data[2]
        dxim = data[4]
        return x, xip, xim, dxip, dxim

    def make(self):
        
        if self.ns:
            print('Processing gammat')
            self.make_gammat()
        if self.nn:
            print('Processing wtheta')
            self.make_wt()
        if self.ss:
            print('Processing xipm')
            self.make_xipm()

        return 0

    def make_xipm(self):

        rows, cols = self.nsrc, self.nsrc
        count = 0

        for i in range(self.nsrc):
            for j in range(self.nsrc):
                count+=1
                if j>i:
                    continue

                print(i,j)
                x, xip,xim, dp, dm = self.get_spectra(j,i)

                #posp = positions[(i+1,j+1,"+")]
                ax = plt.subplot(rows,cols,count)

                ax.annotate("(%d, %d)"%(i+1,j+1), (0.2,4.9), textcoords='data', fontsize=9 )
                ax.yaxis.set_tick_params(which='minor', left='off', right='off')
                plt.ylim(-0.5,6.8)
                #plt.yscale("log")
                plt.xscale("log")
                plt.yticks(fontsize=9)
                plt.xticks([10,100],["10", "100"], fontsize=9)

                if j==self.nsrc-1:
                    #plt.ylabel(r"$\theta \xi_+(\theta)$ / arcmin", fontsize=12)
                    plt.xlabel(r"$\theta$ / arcmin", fontsize=8)
                    plt.yticks(fontsize=9)
                    

                if j==0:
                    plt.yticks([0,2,4,6],['0', '2', '4', '6'])
                    if i==self.nsrc-2:
                        plt.ylabel(r"$\theta \xi_\pm(\theta) / \times 10^{-4}$ arcmin", fontsize=9)
                    if i==self.nsrc-1:
                        plt.xticks([0.1, 10,100],["0.1", "10", "100"], fontsize=9)


                else:
                    plt.yticks(visible=False)

                plt.xlim(0.1,270)
                
                plt.axhline(0, color='k', ls=':')

                plt.errorbar(x, 1e4 * x * xip, yerr=dp, ls="none", marker="D", markersize=2,  ecolor="purple", markeredgecolor="purple", markerfacecolor="none", label=r'$\xi_+$')
                plt.errorbar(x, 1e4 * x * xim, yerr=dm, ls="none", marker="o", markersize=2, ecolor="pink", markeredgecolor="pink", markerfacecolor="pink", label=r'$\xi_-$')

                if (i==0) and (j==0):
                    plt.legend(bbox_to_anchor=(1.8, 1), fontsize=9)


        plt.subplots_adjust(hspace=0,wspace=0.0, bottom=0.14, top=0.97, left=0.14)
        plt.savefig("%s/plots/cornerplot-shear_shear.pdf"%self.basedir)
        plt.close()

        return 0

    def make_gammat(self):

        rows, cols = self.nsrc, self.nlens
        count = 0

        for i in range(self.nsrc):
            for j in range(self.nlens):
                count+=1

                print(i,j)
                x, gt, _, dgt, _ = self.get_spectra(j,i, ptype='position_shear')

                #posp = positions[(i+1,j+1,"+")]
                ax = plt.subplot(rows,cols,count)

                ax.annotate("(%d, %d)"%(i+1,j+1), (0.2,0.6), textcoords='data', fontsize=9 )
                ax.yaxis.set_tick_params(which='minor', left='off', right='off')
                plt.ylim(-0.55,1)
                #plt.yscale("log")
                plt.xscale("log")
                plt.yticks(fontsize=9)
                plt.xticks([10,100],["10", "100"], fontsize=9)

                if j==self.nlens-1:
                    #plt.ylabel(r"$\theta \xi_+(\theta)$ / arcmin", fontsize=12)
                    plt.xlabel(r"$\theta$ / arcmin", fontsize=8)
                    plt.yticks(fontsize=9)

                if j==0:
                    plt.yticks([0, 0.5, 1.])
                    if i==self.nsrc-2:
                        plt.ylabel(r"$\theta \gamma_t(\theta) / \times 10^{-2}$ arcmin", fontsize=9)
                    if i==self.nsrc-1:
                        plt.xticks([0.1, 10,100],["0.1", "10", "100"], fontsize=9)

                else:
                    plt.yticks(visible=False)

                plt.xlim(0.1,270)
                
                plt.axhline(0, color='k', ls=':')

                plt.errorbar(x, 1e2 * x * gt, yerr=dgt, ls="none", marker="D", markersize=2,  ecolor="purple", markeredgecolor="purple", markerfacecolor="none", label=r'$\xi_+$')

                if (i==0) and (j==0):
                    plt.legend(bbox_to_anchor=(1.8, 1), fontsize=9)


        plt.subplots_adjust(hspace=0,wspace=0.0, bottom=0.14, top=0.97, left=0.14)
        plt.savefig("%s/plots/cornerplot-position_shear.pdf"%self.basedir)
        plt.close()

        return 0

    def make_wt(self):

        rows, cols = self.nlens, self.nlens
        count = 0

        for i in range(self.nlens):
            for j in range(self.nlens):
                count+=1
                if j>i:
                    continue

                print(i,j)
                x, wt, _, dwt, _ = self.get_spectra(j,i, ptype='position_position')

                #posp = positions[(i+1,j+1,"+")]
                ax = plt.subplot(rows,cols,count) 

                ax.annotate("(%d, %d)"%(i+1,j+1), (0.2,2), textcoords='data', fontsize=9 )
                ax.yaxis.set_tick_params(which='minor', left='off', right='off')
                plt.ylim(-1,3)
                #plt.yscale("log")
                plt.xscale("log")
                plt.yticks(fontsize=9)
                plt.xticks([10,100],["10", "100"], fontsize=9)

                if j==self.nlens-1:
                    #plt.ylabel(r"$\theta \xi_+(\theta)$ / arcmin", fontsize=12)
                    plt.xlabel(r"$\theta$ / arcmin", fontsize=8)
                    plt.yticks(fontsize=9)

                if j==0:
                    plt.yticks([-1,0,1,2])
                    if i==self.nlens-1:
                        #plt.ylabel(r"$\theta w(\theta)$ / arcmin", fontsize=9)
                        plt.ylabel(r"$\theta w(\theta)$", fontsize=9)
                        plt.xticks([0.1, 10,100],["0.1", "10", "100"], fontsize=9)

                else:
                    plt.yticks(visible=False)

                plt.xlim(0.1,270)   
                plt.axhline(0, color='k', ls=':')

                plt.errorbar(x, x*wt, yerr=dwt, ls="none", marker="D", markersize=2,  ecolor="purple", markeredgecolor="purple", markerfacecolor="none", label=r'$\xi_+$')


        plt.subplots_adjust(hspace=0,wspace=0.0, bottom=0.14, top=0.97, left=0.14)
        plt.savefig("%s/plots/cornerplot-position_position.pdf"%self.basedir)
        plt.close()

        return 0


def main():

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('config', type=str, action='store')
    parser.add_argument('--xipm', '-ss',action='store_true')
    parser.add_argument('--gammat', '-ns',action='store_true')
    parser.add_argument('--wtheta', '-nn',action='store_true')
    args = parser.parse_args()

    config = yaml.load(open(args.config,'rb'))

    plotter = cornerplot(config, args)
    plotter.make()

main()