import os
import FWCore.ParameterSet.Config as cms

# ZDC energy for 2023 PbPb.
#
# 2023 cannot use reco::Centrality: PromptReco carries no centrality object at
# all, and the 14Feb2025 rereco carries one whose ZDC sums show no neutron peak
# structure.  So the sum is taken from the QIE10 digis with the hardcoded 2023
# calibration (TS2-TS1, EM weight 0.1, per-side 0.9397 / 0.5031).
#
# CRITICAL: the 2023 ZDC channel mapping is NOT in the global tag.  Without the
# custom electronics map below, the digis carry wrong HcalZDCDetIds and every
# channel selection -- by DetId or by position -- reads the wrong channels and
# still returns a plausible-looking number.  es_prefer makes the text map win
# over the database one.

process = cms.Process('ZDC2023')
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 20000
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Prompt_v7', '')

# ---------------- 2023 ZDC electronics map ----------------
from CondCore.CondDB.CondDB_cfi import *
# NOTE vs the 2023-era cfg: PoolDBESSource in CMSSW_16_1_x no longer accepts
# 'timetype' or 'authenticationMethod' -- passing them is a fatal Configuration
# error.  Dropped; the defaults are what those values asked for anyway.
process.es_pool = cms.ESSource("PoolDBESSource",
    toGet = cms.VPSet(
        cms.PSet(
            record = cms.string("HcalElectronicsMapRcd"),
            tag = cms.string("HcalElectronicsMap_2021_v2.0_data")
        )
    ),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS')
)
process.es_prefer = cms.ESPrefer('HcalTextCalibrations', 'es_ascii')
process.es_ascii = cms.ESSource(
    'HcalTextCalibrations',
    input = cms.VPSet(
        cms.PSet(
            object = cms.string('ElectronicsMap'),
            file = cms.FileInPath("FlowCorrAna/DiHadronCorrelationAnalyzer/data/emap_2023_newZDC_v3.txt")
        )
    )
)
# ----------------------------------------------------------

_redir = os.environ.get('REDIR','root://cms-xrd-global.cern.ch/')
_flist = os.environ.get('FILELIST','')
FILES = [(_redir + l.strip()) if l.strip().startswith('/store') else l.strip()
         for l in open(_flist)] if _flist else []
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILES))
# Read only what the analyzer needs.  Two reasons: 2023 MINIAOD files differ in
# branch content between runs (some carry slimmedDisplacedMuons collections and
# some do not), which makes PoolSource abort a multi-file job with
# 'MismatchedInputFiles'; and dropping the rest is much faster.
process.source.inputCommands = cms.untracked.vstring(
    'drop *',
    'keep *_hiCentrality_*_*',
    'keep *_hcalDigis_*_*',
    'keep edmTriggerResults_*_*_HLT',
    'keep *_offlineSlimmedPrimaryVertices_*_*',
)
process.source.dropDescendantsOfDroppedBranches = cms.untracked.bool(False)
# Reading over xrootd, an unreadable file otherwise aborts the whole job with a
# FatalRootError before any event is processed.  Skip it and keep going.
process.source.skipBadFiles = cms.untracked.bool(True)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(int(os.environ.get('MAXEV','-1'))))

import FWCore.PythonUtilities.LumiList as LumiList
_json = os.environ.get('GOLDEN_JSON',
    '/afs/cern.ch/user/d/davidlw/json/Cert_Collisions2023HI_374288_375823_Golden.json')
if os.path.exists(_json):
    process.source.lumisToProcess = LumiList.LumiList(filename=_json).getVLuminosityBlockRange()
    print('ZDC_2023: golden JSON %s' % os.path.basename(_json))

# The SAME frozen bit scheme as ZDC_MultiYear_cfg.py, so trigBits means the same
# thing in every year and masks are portable.  Read off the 2023 menu:
#
#   bit 0  ZeroBias_MaxPixelCluster10000               ABSENT in 2023
#   bit 1  ZeroBias_MinPixelCluster400_MaxPixel10000   present
#   bit 2  ZeroBias_SinglePixelTrackLowPt_MaxPix400    present
#   bit 3  ZDC1nOR_MaxPixelCluster10000                ABSENT in 2023
#   bit 4  ZDC1nOR_MinPixelCluster400_MaxPixel10000    present
#   bit 5  ZDC1nOR_SinglePixelTrackLowPt_MaxPix400     present
#   bits 6-8  the three SingleMuOpen paths             present
#
# The two MaxPixelCluster10000 triggers were introduced after 2023.  They are
# the near-unbiased ones that carried the 2024-2026 measurements, so a 2023
# sample is NOT equivalent: both surviving ZB/ZDCOR variants require central
# pixel activity and therefore bias the neutron-multiplicity mixture.  For a
# like-for-like comparison restrict EVERY year to bits 1,2,4,5 (mask 0x036).
DISCOVER = os.environ.get('DISCOVER','0') == '1'
_DEFAULT_PATHS = [
    'HLT_HIUPC_ZeroBias_MaxPixelCluster10000_v*',
    'HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v*',
    'HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*',
    'HLT_HIUPC_ZDC1nOR_MaxPixelCluster10000_v*',
    'HLT_HIUPC_ZDC1nOR_MinPixelCluster400_MaxPixelCluster10000_v*',
    'HLT_HIUPC_ZDC1nOR_SinglePixelTrackLowPt_MaxPixelCluster400_v*',
    'HLT_HIUPC_SingleMuOpen_BptxAND_MaxPixelCluster1000_v*',
    'HLT_HIUPC_SingleMuOpen_NotMBHF2AND_v*',
    'HLT_HIUPC_SingleMuOpen_NotMBHF2AND_MaxPixelCluster1000_v*',
    'HLT_HIUPC_SingleMuOpen_BptxAND_MaxPixelCluster10000_v*',   # 2026 only
]
_env = [p for p in os.environ.get('PATHS','').split(',') if p]
PATHS = _env if _env else _DEFAULT_PATHS

process.load("FlowCorrAna.DiHadronCorrelationAnalyzer.zdcenergyanalyzer_cfi")
process.zdc_ana.VertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
process.zdc_ana.triggerPathNames = cms.vstring([p.rstrip('*') for p in PATHS])
process.zdc_ana.saveTree = cms.untracked.bool(not DISCOVER)
# the 2023 route
process.zdc_ana.doDigiSum  = cms.untracked.bool(True)
process.zdc_ana.zdcDigiSrc = cms.untracked.InputTag('hcalDigis','ZDC')
process.zdc_ana.requireCentrality = cms.untracked.bool(False)

if DISCOVER:
    process.zdc_ana.dumpTriggerPattern = cms.untracked.string(os.environ.get('DUMPPAT','HIUPC'))
    process.ana_step = cms.Path(process.zdc_ana)
else:
    import HLTrigger.HLTfilters.hltHighLevel_cfi
    process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
    process.hltFilter.andOr = cms.bool(True)
    process.hltFilter.throw = cms.bool(False)
    process.hltFilter.HLTPaths = PATHS if PATHS else ['*']
    process.ana_step = cms.Path(process.hltFilter * process.zdc_ana)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(os.environ.get('OUTFILE','ZDC_2023.root')))
process.schedule = cms.Schedule(process.ana_step)
