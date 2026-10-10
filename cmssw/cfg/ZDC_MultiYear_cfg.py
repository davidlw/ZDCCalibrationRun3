import os
import FWCore.ParameterSet.Config as cms

# ZDC energy vs trigger, for 2024 / 2025 / 2026 PbPb with one configuration.
#
# Paths are matched by SUBSTRING with the trailing "_v" kept, so a single list
# works across years even though the version suffixes differ (e.g. the
# ZeroBias_MaxPixelCluster10000 path is _v13 in 2025 and _v6 in 2026).
#
# Bit assignments are FROZEN so that masks already used on 2025 output stay
# valid; anything new appends at the end:
#   bits 0-2  ZB      -> mask 0x007
#   bits 3-5  ZDCOR   -> mask 0x038
#   bits 6-8  MUON    -> mask 0x1C0
#   bit  9    SingleMuOpen_BptxAND_MaxPixelCluster10000 (new in 2026) -> 0x200
#
# Driven by environment: YEAR, FILELIST, OUTFILE, MAXEV, REDIR.

YEAR    = os.environ.get('YEAR', '2025')
JSONDIR = '/afs/cern.ch/user/d/davidlw/json'
JSON = {
  '2023': 'Cert_Collisions2023HI_374288_375823_Golden.json',
  '2024': 'Cert_Collisions2024_HI_387853_388784_Golden.json',
  '2025': 'Cert_Collisions2025_HI_399465_400426_Golden.json',
  '2026': 'Collisions26HI_5p36TeV_404337_404926_golden.json',
}

process = cms.Process('ZDCANA')
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 50000
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# No GlobalTag, geometry or magnetic field: event products only.

_redir = os.environ.get('REDIR', 'root://cms-xrd-global.cern.ch/')
_flist = os.environ.get('FILELIST', '')
FILES = [(_redir + l.strip()) if l.strip().startswith('/store') else l.strip()
         for l in open(_flist)] if _flist else []
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILES))
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(int(os.environ.get('MAXEV', '-1'))))

import FWCore.PythonUtilities.LumiList as LumiList
_json = os.path.join(JSONDIR, JSON.get(YEAR, ''))
if os.path.exists(_json):
    process.source.lumisToProcess = LumiList.LumiList(filename=_json).getVLuminosityBlockRange()
    print('ZDC_MultiYear: %s golden JSON %s' % (YEAR, os.path.basename(_json)))
else:
    print('ZDC_MultiYear: WARNING no JSON for year %s -- processing everything' % YEAR)

ZB_PATHS = [
    'HLT_HIUPC_ZeroBias_MaxPixelCluster10000_v*',
    'HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v*',
    'HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*',
]
ZDCOR_PATHS = [
    'HLT_HIUPC_ZDC1nOR_MaxPixelCluster10000_v*',
    'HLT_HIUPC_ZDC1nOR_MinPixelCluster400_MaxPixelCluster10000_v*',
    'HLT_HIUPC_ZDC1nOR_SinglePixelTrackLowPt_MaxPixelCluster400_v*',
]
MUON_PATHS = [
    'HLT_HIUPC_SingleMuOpen_BptxAND_MaxPixelCluster1000_v*',
    'HLT_HIUPC_SingleMuOpen_NotMBHF2AND_v*',
    'HLT_HIUPC_SingleMuOpen_NotMBHF2AND_MaxPixelCluster1000_v*',
]
NEW_PATHS = [
    'HLT_HIUPC_SingleMuOpen_BptxAND_MaxPixelCluster10000_v*',   # 2026, absent earlier
]
ALL_PATHS = ZB_PATHS + ZDCOR_PATHS + MUON_PATHS + NEW_PATHS

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)     # tolerate paths absent in a given year
process.hltFilter.HLTPaths = ALL_PATHS

process.load("FlowCorrAna.DiHadronCorrelationAnalyzer.zdcenergyanalyzer_cfi")
process.zdc_ana.VertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
process.zdc_ana.triggerPathNames = cms.vstring([p.rstrip('*') for p in ALL_PATHS])
process.zdc_ana.saveTree = cms.untracked.bool(True)

# ---- ZDC source: centrality for 2024-2026, hardcoded digis for 2023 ONLY ----
#
# 2024-2026 use reco::Centrality.  Its zdcSumPlus/Minus are the calibrated Run 3
# reconstruction (template fit over bxTs={0,2,4} with conditions pedestals), and
# are identical to the forest ZDCRecHitAnalyzerHC sum -- the RPD and EM-channel>5
# cuts are applied upstream in ZdcHitReconstructor_Run3.
#
# 2023 CANNOT use it: PromptReco carries no centrality object at all, and the
# 14Feb2025 rereco carries one whose sums show no neutron peak structure.  So
# 2023 is summed straight from the QIE10 digis with the hardcoded 2023
# calibration (TS2-TS1, EM weight 0.1, per-side 0.9397/0.5031).
#
# The digi sum is a CRUDER estimator even when correct -- validated against
# centrality on 2025 it correlates only 0.88-0.90, because TS2-TS1 has no
# conditions pedestal and no out-of-time subtraction.  So it is a fallback for
# 2023, never a substitute where centrality is available.
#
# Set DODIGI=1 to additionally record it on another year purely for validation.
DIGI_ONLY_YEAR = '2023'
_wantDigi = (YEAR == DIGI_ONLY_YEAR) or (os.environ.get('DODIGI','0') == '1')
if _wantDigi:
    process.zdc_ana.doDigiSum  = cms.untracked.bool(True)
    process.zdc_ana.zdcDigiSrc = cms.untracked.InputTag('hcalDigis','ZDC')
    # 2023 has no usable centrality object, so do not drop those events
    process.zdc_ana.requireCentrality = cms.untracked.bool(YEAR != DIGI_ONLY_YEAR)
    print('ZDC_MultiYear: %s digi sum ON (%s)' %
          (YEAR, 'primary, 2023 route' if YEAR == DIGI_ONLY_YEAR else 'validation alongside centrality'))
if YEAR != DIGI_ONLY_YEAR:
    print('ZDC_MultiYear: %s ZDC from reco::Centrality' % YEAR)

process.ana_step = cms.Path(process.hltFilter * process.zdc_ana)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(os.environ.get('OUTFILE', 'ZDC_%s.root' % YEAR)))
process.schedule = cms.Schedule(process.ana_step)
