import FWCore.ParameterSet.Config as cms

zdc_ana = cms.EDAnalyzer('ZDCEnergyAnalyzer',

  # ZDC comes from reco::Centrality (zdcSumPlus/zdcSumMinus), which is the
  # calibrated sum for 2024-2026 data.  2023 is the special case that needs the
  # emap + zdcdigi + QWzdcreco chain; nothing here does.
  centralitySrc = cms.InputTag("hiCentrality"),

  # Only read when useCentralityBin is True.  Not needed for the ZDC spectra.
  centralityBinLabel = cms.InputTag("centralityBin", "HFtowers"),

  TriggerResultCollection = cms.InputTag("TriggerResults", "", "HLT"),

  # Recorded, never cut on.  On MINIAOD point this at the slimmed collection.
  VertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),

  # The paths whose decision is stored.  Bit i of trigBits in the tree is entry
  # i of this list, and the per-trigger histograms follow the same order, so
  # keep the order stable once you have started producing files.
  #
  # Defaults are the 2025 PbPb UPC paths.  The cfg overrides this from the same
  # list it gives hltHighLevel, so the filter and the bookkeeping cannot drift
  # apart; edit the cfg, not this, when the menu changes.  Verify against a new
  # period's menu before trusting output: hltHighLevel with throw=False accepts
  # nothing when no configured path exists, which looks exactly like a real
  # empty spectrum.  Run the cfg with DISCOVER_TRIGGERS=True to print the menu.
  triggerPathNames = cms.vstring(
      # zero bias                                        -> trigBits 0-2
      'HLT_HIUPC_ZeroBias_MaxPixelCluster10000_v',
      'HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v',
      'HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v',
      # ZDC OR                                           -> trigBits 3-5
      'HLT_HIUPC_ZDC1nOR_MaxPixelCluster10000_v',
      'HLT_HIUPC_ZDC1nOR_MinPixelCluster400_MaxPixelCluster10000_v',
      'HLT_HIUPC_ZDC1nOR_SinglePixelTrackLowPt_MaxPixelCluster400_v',
      # single muon                                      -> trigBits 6-8
      'HLT_HIUPC_SingleMuOpen_BptxAND_MaxPixelCluster1000_v',
      'HLT_HIUPC_SingleMuOpen_NotMBHF2AND_v',
      'HLT_HIUPC_SingleMuOpen_NotMBHF2AND_MaxPixelCluster1000_v',
  ),

  # The tree is what the analysis is done from; the histograms are there so a
  # job can be checked without an offline step.  Turn the tree off only once
  # the binning is settled and you are running over everything.
  saveTree = cms.untracked.bool(True),

  saveVertexInfo = cms.untracked.bool(True),
  # Slimmed vertices carry no track references, so tracksSize() is 0 for every
  # vertex in MINIAOD.  Leave at 0 there; use 2 on AOD to match SoS's counting.
  minVtxTracks = cms.untracked.int32(0),

  useCentralityBin = cms.untracked.bool(False),

  # Monitoring-histogram ranges only -- the tree keeps the full value, so a bad
  # guess here costs nothing.
  #
  # Narrow range, 1 unit/bin: this is where the events are.  Measured on run
  # 400387, the ZDC+ distribution has median 41 and 75th percentile 97, so a
  # range of 0-100000 (the first guess) put the entire 0n/1n region into its
  # first four bins.  Starts below zero so pedestal-negative values show up
  # instead of piling into the underflow bin.
  zdcHistBins = cms.untracked.int32(2200),
  zdcHistMin  = cms.untracked.double(-200.),
  zdcHistMax  = cms.untracked.double(2000.),
  # Bins per axis of every 2D.  Cost is bins^2, and with one 2D per trigger
  # these dominate a histogram-only output file: 500 -> 300 cuts them by 2.8x.
  zdc2DBins = cms.untracked.int32(300),
  # Per-trigger 2D correlation plots.  Turn off for a large multi-PD run if the
  # inclusive 2D is enough; the per-trigger 1D spectra are kept either way.
  saveTrigger2D = cms.untracked.bool(True),
  # All fills have unit weight, so GetBinError() is sqrt(N) with or without the
  # sumw2 array -- storing it just doubles the file.  Turn ON if a weighted
  # fill is ever added, otherwise the errors would be wrong.
  storeErrors = cms.untracked.bool(False),

  # Tree slimming (see the header).  The event number and hfSumET were never
  # read by any macro and cost 43% of the tree; ZDC sums rounded to integers
  # compress 3x better and every reader still sees a Float_t.  nPixel and bx
  # are 16-bit.  compression is alg*100+level: 505 = ZSTD level 5,
  # 101 = the old zlib default.  Net: 15.5 -> 4.8 bytes/event.
  storeEvent  = cms.untracked.bool(False),
  storeHF     = cms.untracked.bool(False),
  roundZDC    = cms.untracked.bool(True),
  compression = cms.untracked.int32(505),

  # Wide range for the Xn tail (99th percentile ~45000, seen up to ~103000).
  # Inclusive histograms only; per-trigger tails come from the tree.
  zdcWideBins = cms.untracked.int32(2000),
  zdcWideMax  = cms.untracked.double(200000.),

  # Non-empty: label-fill every fired path containing this substring into
  # hAllTrig, so a short unfiltered run prints the actual menu.
  dumpTriggerPattern = cms.untracked.string(''),
)
