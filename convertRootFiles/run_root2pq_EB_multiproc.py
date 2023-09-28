import os, glob, re
import shutil
from multiprocessing import Pool

def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [int(c) if c.isdigit() else c for c in re.split('([0-9]+)',s)]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)

def run_process(process):
    os.system('python %s'%process)

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-d', '--genDR', default=10, type=int, help='gen-level dR.')
parser.add_argument('-p', '--p_drop', default=1.00, type=float, help='p(drop) scale.')
#parser.add_argument('-w', '--wgt_file', default=None, type=str, help='Weight file.')
parser.add_argument('-w', '--wgt_files', default=None, nargs='+', type=str, help='Weight file.')
parser.add_argument('-b', '--batch_size', default=1, type=int, help='Batch size.')
args = parser.parse_args()

genDR = args.genDR
p_drop_scale = args.p_drop
wgt_files = args.wgt_files

#xrootd='root://cmsxrootd.fnal.gov' # FNAL
xrootd='root://cmseos.fnal.gov' # FNAL
#xrootd='root://eoscms.cern.ch' # CERN
eosDir='/eos/uscms/store/user/lpcml/mandrews/IMG'

decay = 'GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8_MINIAODSIM' # all
decay = 'GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_MINIAODSIM' # 100
decay = 'QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_MINIAODSIM' # 163
decay = 'DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8_MINIAODSIM' # 163
date_str = '190912_201036'
date_str = '190912_200744'
date_str = '190912_200901'
date_str = '190912_210041'
decay = 'h24gamma_1j_1M_100MeV_PU2016_MINIAODSIM'
decay = 'h24gamma_1j_1M_200MeV_PU2016_MINIAODSIM'
#decay = 'h24gamma_1j_1M_400MeV_PU2016_MINIAODSIM'
date_str = '190922_185712'
date_str = '190922_221716'
#date_str = '190922_185411'

# Paths to input files
#rhFileList = '%s/%s/2016_%s/%s/*/output_*.root'%(eosDir, decay.partition('_MINIAODSIM')[0], decay, date_str)
rhFileList = '%s/h24gamma_2016/%s/%s/*/output_*.root'%(eosDir, decay, date_str)
print(" >> Input file list: %s"%rhFileList)
rhFileList = glob.glob(rhFileList)
assert len(rhFileList) > 0
print(" >> %d files found"%len(rhFileList))
rhFileList = [('%s/%s'%(xrootd, rhFile)).replace('/eos/uscms','') for rhFile in rhFileList]
print(' >> Input File[0]: %s'%rhFileList[0])
sort_nicely(rhFileList)

# Weights file
if wgt_files is not None:
  #wgt_file = '%s_mvpt_weights_pdrop%.2f.npz'%(decay, p_drop_scale)
  #wgt_files = args.wgt_files
  for wgt_file in wgt_files:
      assert os.path.isfile(wgt_file)

# Output path
outDir='/uscms/physics_grp/lpcml/nobackup/mandrews' # NOTE: Space here is limited, transfer files to EOS after processing
outDir='%s/%s'%(outDir, decay)
if not os.path.isdir(outDir):
    os.makedirs(outDir)
print(' >> Output directory: %s'%outDir)

proc_file = 'convert_root2pq_EB.py'
if wgt_files is not None:
    #processes = ['%s -i %s -o %s -d %s -n %d -w %s'%(proc_file, rhFile, outDir, decay, i+1, wgt_files) for i,rhFile in enumerate(rhFileList)]
    processes = ['%s -i %s -o %s -d %s -n %d -w %s'%(proc_file, rhFile, outDir, decay, i+1, ' '.join(wgt_files)) for i,rhFile in enumerate(rhFileList)]
else:
  #processes = ['%s -i %s -o %s -d %s -n %d'%(proc_file, rhFile, outDir, decay, i+1) for i,rhFile in enumerate(rhFileList)]
  #processes = ['%s -i %s -o %s -d %s'%(proc_file, ' '.join(rhFileList), outDir, decay)]
  processes = []
  for it,i in enumerate(range(0, len(rhFileList), args.batch_size)):
    rhFileList_batch = rhFileList[i:i+args.batch_size]
    #print(it, i, i+args.batch_size, len(rhFileList_batch))
    processes.append('%s -i %s -o %s -d %s -n %d'%(proc_file, ' '.join(rhFileList_batch), outDir, decay, it))
#print(' >> Process[0]: %s'%processes[0])

#os.system('python %s -i %s -o %s -d %s -n %d -w %s'%(proc_file, rhFileList[0], outDir, decay, 1, wgt_file))
pool = Pool(processes=len(processes))
pool.map(run_process, processes)
pool.close()
pool.join()

for f in glob.glob('%s/%s*.parquet.*'%(outDir, decay)):
  pass
  #shutil.move(f, '%s/%s/2016_%s/%s/'%(eosDir, decay.partition('_MINIAODSIM')[0], decay, date_str))
  shutil.move(f, '%s/h24gamma_2016/%s/%s/'%(eosDir, decay, date_str))
