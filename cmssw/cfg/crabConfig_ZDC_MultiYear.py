from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from http.client import HTTPException
from CRABClient.UserUtilities import config
import os, sys

# Slim ZDC trees, one full PD per year, 2023-2026 PbPb MINIAOD.
#
# The PSet (ZDC_MultiYear_cfg.py / ZDC_2023_MINIAOD_cfg.py) is driven by
# environment variables.  CRAB evaluates the PSet HERE, at submission, and
# ships the pickled process, so setting os.environ before each submit is all
# that is needed -- nothing is read from the environment on the worker.
# FILELIST is left empty so CRAB fills the input files itself.
#
# 2024 is THREE datasets covering disjoint run ranges (A-v1 387757-388192,
# B-v1 388298-388459, B-v2 388468-388784).  They are submitted as three tasks
# and merged offline; the ZDC calibration eras do NOT follow the dataset
# boundaries (see ZDC_STUDY_SUMMARY.md sec. 20), so the era split is done on
# run number after the merge.
#
# 2026 routes ZDC1nOR only to HIForward22, so that PD is needed as well.
#
# Run:   source /cvmfs/cms.cern.ch/common/crab-setup.sh
#        voms-proxy-init -voms cms -valid 192:00 ; (myproxy delegation prompt)
#        python3 crabConfig_ZDC_MultiYear.py [task ...]
# with no arguments every task in TASKS is submitted; with arguments only the
# named ones.  Each task is submitted in its OWN python process (see the end):
# a second crabCommand('submit') in one process re-imports the same PSet
# module and dies silently after creating the project directory.

DATE_TAG = '20260915v1'
JSONDIR  = '/afs/cern.ch/user/d/davidlw/json'
JSON = {
  '2023': 'Cert_Collisions2023HI_374288_375823_Golden.json',
  '2024': 'Cert_Collisions2024_HI_387853_388784_Golden.json',
  '2025': 'Cert_Collisions2025_HI_399465_400426_Golden.json',
  '2026': 'Collisions26HI_5p36TeV_404337_404926_golden.json',
}

# name, dataset, year, pset
TASKS = [
  # 2023: ONLY PromptReco-v2 (runs 374668-375831, 78 of the 93 golden runs) is
  # used.  PromptReco-v1 (373870-374667, 15 golden runs at the low-luminosity
  # start of the run) is not part of the analysis.
  ('ZDC_2023v2_HIForward1', '/HIForward1/HIRun2023A-PromptReco-v2/MINIAOD',  '2023', 'ZDC_2023_MINIAOD_cfg.py'),
  ('ZDC_2024Av1_HIForward1','/HIForward1/HIRun2024A-PromptReco-v1/MINIAOD',  '2024', 'ZDC_MultiYear_cfg.py'),
  ('ZDC_2024Bv1_HIForward1','/HIForward1/HIRun2024B-PromptReco-v1/MINIAOD',  '2024', 'ZDC_MultiYear_cfg.py'),
  ('ZDC_2024Bv2_HIForward1','/HIForward1/HIRun2024B-PromptReco-v2/MINIAOD',  '2024', 'ZDC_MultiYear_cfg.py'),
  ('ZDC_2025_HIForward1',   '/HIForward1/HIRun2025A-PromptReco-v1/MINIAOD',  '2025', 'ZDC_MultiYear_cfg.py'),
  # 2026 needs ONE PD.  An early 40k-event sample suggested ZDC1nOR was routed
  # only to HIForward22, but that sample landed on commissioning runs outside
  # the certified range.  Full processing shows HIForward1 and HIForward22 carry
  # the same menu at the same rates (1.06G ZDC1nOR events each), so HIForward22
  # is an independent cross-check sample, not a requirement.
  ('ZDC_2026_HIForward1',   '/HIForward1/HIRun2026A-PromptReco-v1/MINIAOD',  '2026', 'ZDC_MultiYear_cfg.py'),
]

# ~13k events per LS per PD, 50 LS/job -> ~650k events, ~5 MB slim output, ~6 min.
UNITS_PER_JOB = 50

config = config()
config.section_("General")
config.General.workArea = '/eos/user/d/davidlw/crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB = 2500

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = UNITS_PER_JOB
config.Data.publication = False
config.Data.allowNonValidInputDataset = True
config.Data.outLFNDirBase = '/store/group/phys_heavyions/davidlw/ZDCslim/'

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

def submit(config):
    try:
        crabCommand('submit', config = config, dryrun=False)
    except HTTPException as hte:
        print("Failed submitting task: %s" % (hte.headers))
    except ClientException as cle:
        print("Failed submitting task: %s" % (cle))

wanted = set(sys.argv[1:])
if len(wanted) != 1:
    # fan out: one subprocess per task
    import subprocess
    for name, dataset, year, pset in TASKS:
        if wanted and name not in wanted: continue
        subprocess.call([sys.executable, __file__, name])
    sys.exit(0)

for name, dataset, year, pset in TASKS:
    if wanted and name not in wanted: continue
    json = os.path.join(JSONDIR, JSON[year])
    # PSet environment, read when CRAB evaluates the PSet at submission
    os.environ['YEAR'] = year
    os.environ['FILELIST'] = ''
    os.environ['MAXEV'] = '-1'
    os.environ['OUTFILE'] = 'ZDC_%s.root' % year
    os.environ['GOLDEN_JSON'] = json        # used by the 2023 cfg
    config.General.requestName = '%s_%s' % (name, DATE_TAG)
    config.JobType.psetName = pset
    config.JobType.inputFiles = [json]
    config.JobType.outputFiles = ['ZDC_%s.root' % year]
    config.Data.inputDataset = dataset
    config.Data.lumiMask = json
    config.Data.outputDatasetTag = config.General.requestName
    print("Submitting %s : %s  (year %s, %s)" % (name, dataset, year, pset))
    submit(config)
