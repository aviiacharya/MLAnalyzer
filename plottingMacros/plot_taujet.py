from torch.utils.data import *
from sklearn.metrics import roc_curve, auc

import pyarrow.parquet as pq
import pyarrow as pa # pip install pyarrow==0.7.1
import ROOT
import numpy as np
np.random.seed(0)
import glob, os

import dask.array as da

from scipy.misc import imresize

import matplotlib.pyplot as plt
#%matplotlib inline
from matplotlib.colors import LogNorm, ListedColormap, LinearSegmentedColormap
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

from skimage.measure import block_reduce
from numpy.lib.stride_tricks import as_strided

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-i', '--infile', default='output_DoubleTau.root', type=str, help='Input root file.')
parser.add_argument('-o', '--outdir', default='HToAAToTauTau', type=str, help='Output directory.')
parser.add_argument('-n', '--nEvents', default=10, type=int, help='Number of events.')
parser.add_argument('-s', '--skipEvents', default=-1, type=int, help='Number of events.')
args = parser.parse_args()

fileStr = args.infile
outDir  = args.outdir
f0s = glob.glob(fileStr)

class ParquetDataset(Dataset):
    def __init__(self, filename):
        self.parquet = pq.ParquetFile(filename)
        self.cols = None # read all columns
        #self.cols = ['X_jet.list.item.list.item.list.item','y'] 
    def __getitem__(self, index):
        data = self.parquet.read_row_group(index, columns=self.cols).to_pydict()
        data['X_jet'] = np.float32(data['X_jet'][0])
        data['y'] = np.float32(data['y'])
        data['m0'] = np.float32(data['m0'])
        data['pt'] = np.float32(data['pt'])
        # Preprocessing
        data['X_jet'][data['X_jet'] < 1.e-3] = 0. # Zero-Suppression
        data['X_jet'][-5,...] = 25.*data['X_jet'][-5,...] # For HCAL: to match pixel intensity dist of other layers
        data['X_jet'] = data['X_jet']/100. # To standardize
        return dict(data)
    def __len__(self):
        return self.parquet.num_row_groups

def custom_div_cmap(numcolors=11, name='custom_div_cmap',mincol='blue', midcol='white', maxcol='red'):
    cmap = LinearSegmentedColormap.from_list(name=name,colors=[mincol, midcol, maxcol],N=numcolors)
    return cmap

pink_map = custom_div_cmap(50, mincol='#FFFFFF', midcol='#F699CD' ,maxcol='#FF1694')

def plotJet(img, mins, maxs, str_):
    im = plt.imshow(np.zeros_like(img[6,:,:]), cmap='Greys', vmin=0., vmax=1., alpha=0.9)
    #if maxs[-1] > 0 : plt.imshow(img[6,:,:], cmap='Greens', norm=LogNorm(), alpha=0.9, vmin=mins[-1], vmax=maxs[-1])
    if maxs[-2] > 0 : plt.imshow(img[5,:,:], cmap='Purples', norm=LogNorm(), alpha=0.9, vmin=mins[-2], vmax=maxs[-2])
    if maxs[-3] > 0 : plt.imshow(img[4,:,:], cmap='Blues', norm=LogNorm(), alpha=0.9, vmin=mins[-3], vmax=maxs[-3])
    if maxs[-4] > 0 : plt.imshow(img[3,:,:], cmap='Greens', norm=LogNorm(), alpha=0.9, vmin=mins[-4], vmax=maxs[-4])
    if maxs[-5] > 0 : plt.imshow(img[2,:,:], cmap='Greys',  norm=LogNorm(), alpha=0.9, vmin=mins[-5], vmax=maxs[-5])
    if maxs[-6] > 0 : plt.imshow(img[1,:,:], cmap='Blues',  norm=LogNorm(), alpha=0.9, vmin=mins[-6], vmax=maxs[-6])
    if maxs[-7] > 0 : plt.imshow(img[0,:,:], cmap='Oranges',norm=LogNorm(), alpha=0.9, vmin=mins[-7], vmax=maxs[-7])
    #plt.colorbar(fraction=0.046, pad=0.04)

    #X AXIS
    ax = plt.axes()
    plt.xlim([0., 125.+0.])
    plt.xticks(np.arange(0,150,25))
    ax_range_x = np.arange(0,125+25,25)
    ax.set_xticks(ax_range_x)
    ax.set_xticklabels(ax_range_x)
    plt.xlabel(r"$\mathrm{i\varphi}'$", size=28) #28, 30
    ax.xaxis.set_tick_params(direction='in', which='major', length=6.)
    ax.xaxis.set_tick_params(direction='in', which='minor', length=3.)

    #Y AXIS
    plt.ylim([125.+0.,0.])
    plt.yticks(np.arange(150,0,25))
    plt.ylabel(r"$\mathrm{i\eta}'$", size=28) #28, 30
    ax_range_y = np.arange(0,125+25,25)
    ax.set_yticks(ax_range_y)
    ax.set_yticklabels(ax_range_y)
    ax.yaxis.set_tick_params(direction='in', which='major', length=6.)
    ax.yaxis.set_tick_params(direction='in', which='minor', length=3.)

    #LEGEND
    #colors = {1:'tab:orange',2:'tab:blue',3:'tab:grey',4:'tab:green',5:'tab:blue',6:'tab:purple',7:'tab:green'}
    colors = {1:'orange',2:'lightblue',3:'grey',4:'green',5:'blue',6:'purple'}
    #labels = {1:'Track pT',2:'ECAL',3:'HCAL',4:'PXB1',5:'PXB2',6:'PXB3',7:'PXB4'}
    labels = {1:'Track pT',2:'ECAL',3:'HCAL',4:'BPix L1',5:'BPix L2',6:'BPix L3'}
    patches =[mpatches.Patch(color=colors[i],label=labels[i]) for i in colors]
    #plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
    plt.legend(handles=patches, loc='best')
    #plt.savefig(str_, bbox_inches='tight')
    plt.savefig(str_, bbox_inches='tight', format='png')
    #plt.show()
    plt.clf()

def plotJet_PBX12(img, mins, maxs, str_):
    plt.imshow(np.zeros_like(img[4,:,:]), cmap='Greys', vmin=0., vmax=1., alpha=0.9)
    if maxs[-3] > 0 : plt.imshow(img[4,:,:], cmap=pink_map, norm=LogNorm(), alpha=0.9, vmin=mins[-3], vmax=maxs[-3])
    if maxs[-4] > 0 : plt.imshow(img[3,:,:], cmap='Reds', norm=LogNorm(), alpha=0.9, vmin=mins[-4], vmax=maxs[-4])
    #plt.colorbar(fraction=0.046, pad=0.04)
    ax = plt.axes()
    plt.xlim([0., 125.+0.])
    plt.xticks(np.arange(0,150,25))
    plt.xlabel(r"$\mathrm{i\varphi}'$", size=28) #28, 30
    ax.xaxis.set_tick_params(direction='in', which='major', length=6.)
    plt.ylim([0., 125.+0.])
    plt.yticks(np.arange(0,150,25))
    plt.ylabel(r"$\mathrm{i\eta}'$", size=28) #28, 30
    ax.yaxis.set_tick_params(direction='in', which='major', length=6.)
    #plt.savefig(str_, bbox_inches='tight')
    plt.savefig(str_, bbox_inches='tight', format='png')
    plt.clf()
    #plt.show()

def plotJet_PBX(img, mins, maxs, str_):
    plt.imshow(np.zeros_like(img[6,:,:]), cmap='Greys', vmin=0., vmax=1., alpha=0.9)
    if maxs[-1] > 0 : plt.imshow(img[6,:,:], cmap='Greens', norm=LogNorm(), alpha=0.9, vmin=mins[-1], vmax=maxs[-1])
    if maxs[-2] > 0 : plt.imshow(img[5,:,:], cmap='Purples', norm=LogNorm(), alpha=0.9, vmin=mins[-2], vmax=maxs[-2])
    if maxs[-3] > 0 : plt.imshow(img[4,:,:], cmap=pink_map, norm=LogNorm(), alpha=0.9, vmin=mins[-3], vmax=maxs[-3])
    if maxs[-4] > 0 : plt.imshow(img[3,:,:], cmap='Reds', norm=LogNorm(), alpha=0.9, vmin=mins[-4], vmax=maxs[-4])
    if maxs[-7] > 0 : plt.imshow(img[0,:,:], cmap='Oranges',norm=LogNorm(), alpha=0.9, vmin=mins[-7], vmax=maxs[-7])
    #plt.colorbar(fraction=0.046, pad=0.04)
    ax = plt.axes()
    plt.xlim([0., 125.+0.])
    plt.xticks(np.arange(0,150,25))
    plt.xlabel(r"$\mathrm{i\varphi}'$", size=28) #28, 30
    ax.xaxis.set_tick_params(direction='in', which='major', length=6.)
    plt.ylim([0., 125.+0.])
    plt.yticks(np.arange(0,150,25))
    plt.ylabel(r"$\mathrm{i\eta}'$", size=28) #28, 30
    ax.yaxis.set_tick_params(direction='in', which='major', length=6.)
    #plt.savefig(str_, bbox_inches='tight')
    plt.savefig(str_, bbox_inches='tight', format='png')
    plt.clf()
    #plt.show()


def plotJet_chnl(img, cmap_, xmin, xmax, str_):
    plt.imshow(np.zeros_like(img), cmap='Greys', vmin=0., vmax=1., alpha=0.9)
    plt.imshow(img, cmap=cmap_, norm=LogNorm(), alpha=0.9, vmin=xmin, vmax=xmax)
    ax = plt.axes()
    plt.xlim([0., 125.+0.])
    plt.xticks(np.arange(0,150,25))
    plt.xlabel(r"$\mathrm{i\varphi}'$", size=28) #28, 30
    ax.xaxis.set_tick_params(direction='in', which='major', length=6.)
    plt.ylim([0., 125.+0.])
    plt.yticks(np.arange(0,150,25))
    plt.ylabel(r"$\mathrm{i\eta}'$", size=28) #28, 30
    ax.yaxis.set_tick_params(direction='in', which='major', length=6.)
    plt.savefig(str_, bbox_inches='tight', format='png')
    plt.clf()
    #plt.show()

dset_train = ParquetDataset(fileStr)
train_cut = 50
idxs = np.random.permutation(len(dset_train))
train_sampler = sampler.SubsetRandomSampler(idxs[:train_cut])
#train_loader = DataLoader(dataset=dset_train, batch_size=32, num_workers=0, sampler=train_sampler, pin_memory=True)
train_loader = DataLoader(dataset=dset_train, batch_size=2, num_workers=0, shuffle=False, pin_memory=True)
for i, data in enumerate(train_loader):
    if i < args.skipEvents: continue
    print " Event ", i
      
    if i == args.nEvents: break
    X_train = data['X_jet']
    y_train = data['y']

    plt.rcParams["font.family"] = "Helvetica"
    plt.rcParams["figure.figsize"] = (12,12)
    #plt.rcParams["axes.facecolor"] = "white"
    plt.rcParams.update({'font.size': 26})
    
    cmap = ['Oranges','Blues','Greys','Reds',pink_map,'Purples','Greens']
    min_ = 0.0001

    for jet in range(2):
        img = X_train[jet,:,:,:]
        print("JET LABEL IS  ", y_train[jet])
        #Selecting only taus
        if y_train[jet] == 0: continue

        for ch in range(7):
            img_ = img[ch,:,:]
            max_ = img_.max()
            if max_ == 0: continue
            print "Channel ", ch, " , Max = ", max_
            plotJet_chnl(img_, cmap[ch], min_, max_, 'images/%s/tau_event%d_jet%d_chnl%d.png'%(outDir,i,jet,ch))

        mins = [0.0001]*7
        maxs = [X_train[jet,0,:,:].max(), X_train[jet,1,:,:].max(), X_train[jet,2,:,:].max(), X_train[jet,3,:,:].max(), X_train[jet,4,:,:].max(), X_train[jet,5,:,:].max(), X_train[jet,6,:,:].max()]
        print "Min = ", mins, " | Max = ", maxs
        plotJet(img, mins, maxs, 'images/%s/tau_event%d_jet%d.png'%(outDir,i,jet))
        #plotJet_PBX12(img, mins, maxs, 'images/%s/tau_event%d_jet%d_PBX12.png'%(outDir,i,jet))
        #plotJet_PBX(img, mins, maxs, 'images/%s/tau_event%d_jet%d_PBX.png'%(outDir,i,jet))
