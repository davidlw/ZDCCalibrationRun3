#include "FlowCorrAna/DiHadronCorrelationAnalyzer/interface/ZDCEnergyAnalyzer.h"
#include "FlowCorrAna/DiHadronCorrelationAnalyzer/interface/QWZDC2018Helper.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TString.h>

#include <iostream>
#include <cmath>

namespace {
  // trigBits_ is one 64-bit word, so the configured list cannot be longer.
  const size_t kMaxTrig = 64;
}

ZDCEnergyAnalyzer::ZDCEnergyAnalyzer(const edm::ParameterSet& iConfig) :
  hNEvt(0), hTrigCount(0), hAllTrig(0),
  hZdcPlus(0), hZdcMinus(0), hZdcSum(0), hZdcMinusVsPlus(0),
  hZdcPlusWide(0), hZdcMinusWide(0), hZdcSumWide(0), hZdcMinusVsPlusWide(0),
  hHF(0), hNVtx(0),
  zdcTree(0),
  runNb_(0), lsNb_(0), eventNb_(0),
  zdcPlus_(-1.f), zdcMinus_(-1.f), hfSumET_(-1.f),
  nPixel_(65535), bx_(0), nVtx_(-1), centBin_(-1), trigBits_(0),
  nProcessed_(0), nNoCentrality_(0)
{
  // Every fill here has unit weight, so the sumw2 array is exactly the bin
  // content and GetBinError() returns sqrt(N) without it -- storing it doubles
  // the file for no information.  It stays available for anyone who later adds
  // a weighted fill, in which case this MUST be turned back on.
  storeErrors_ = iConfig.getUntrackedParameter<bool>("storeErrors", false);
  if(storeErrors_) TH1::SetDefaultSumw2();
  usesResource("TFileService");

  // Tree slimming.  ZSTD level 5 (505) is ~15% smaller than zlib-1 on this
  // tree at negligible CPU cost; set on the file before any tree is created so
  // the tree inherits it.
  storeEvent_  = iConfig.getUntrackedParameter<bool>("storeEvent", false);
  storeHF_     = iConfig.getUntrackedParameter<bool>("storeHF", false);
  roundZDC_    = iConfig.getUntrackedParameter<bool>("roundZDC", true);
  compression_ = iConfig.getUntrackedParameter<int>("compression", 505);
  theOutputs->file().SetCompressionSettings(compression_);

  token_centrality = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));

  // Hardcoded digi-level ZDC sum.  Constants default to the 2023 values from
  // SoSMeanPtAnalyzer; they are configurable because they are a calibration,
  // not a property of the code.
  doDigiSum_        = iConfig.getUntrackedParameter<bool>("doDigiSum", false);
  requireCentrality_= iConfig.getUntrackedParameter<bool>("requireCentrality", true);
  emWeight_         = iConfig.getUntrackedParameter<double>("emWeight", 0.1);
  calibPlus_        = iConfig.getUntrackedParameter<double>("calibPlus", 0.9397);
  calibMinus_       = iConfig.getUntrackedParameter<double>("calibMinus", 0.5031);
  tsSig_            = iConfig.getUntrackedParameter<int>("tsSignal", 2);
  tsPed_            = iConfig.getUntrackedParameter<int>("tsPedestal", 1);
  if(doDigiSum_)
    token_zdcdigi = consumes<QIE10DigiCollection>(
        iConfig.getUntrackedParameter<edm::InputTag>("zdcDigiSrc", edm::InputTag("hcalDigis","ZDC")));
  token_triggerResults = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResultCollection"));

  triggerNames_ = iConfig.getParameter<std::vector<std::string> >("triggerPathNames");
  if(triggerNames_.size() > kMaxTrig)
    throw cms::Exception("Configuration")
      << "triggerPathNames has " << triggerNames_.size()
      << " entries; the trigger bit word holds at most " << kMaxTrig << ".";

  // The centrality *bin* is a separate product from the centrality object and
  // is not needed for the ZDC spectra, so it is off by default: asking for it
  // on a file that does not carry it would abort the job for nothing.
  useCentralityBin_ = iConfig.getUntrackedParameter<bool>("useCentralityBin", false);
  if(useCentralityBin_)
    token_centralitybin = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinLabel"));

  // Vertex multiplicity is recorded, never cut on -- see the header.
  saveVertexInfo_ = iConfig.getUntrackedParameter<bool>("saveVertexInfo", true);
  if(saveVertexInfo_)
    token_vertices = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("VertexCollection"));
  // Slimmed vertices in MINIAOD carry no track references, so tracksSize() is 0
  // there for every vertex; leave this at 0 unless running on AOD.
  minVtxTracks_ = iConfig.getUntrackedParameter<int>("minVtxTracks", 0);

  saveTree_ = iConfig.getUntrackedParameter<bool>("saveTree", true);
  dumpTriggerPattern_ = iConfig.getUntrackedParameter<std::string>("dumpTriggerPattern", "");

  zdcHistBins_ = iConfig.getUntrackedParameter<int>("zdcHistBins", 3040);
  zdcHistMin_  = iConfig.getUntrackedParameter<double>("zdcHistMin", -200.);
  // Must stay above the neutron peaks: 1n ~2450, 2n ~4800, 3n ~6600.
  zdcHistMax_  = iConfig.getUntrackedParameter<double>("zdcHistMax", 15000.);
  zdcWideBins_ = iConfig.getUntrackedParameter<int>("zdcWideBins", 2000);
  // Each 2D costs bins^2 doubles.  At 500 that is 250k bins per histogram and,
  // with one per trigger, the 2Ds are most of a histogram-only output file.
  zdc2DBins_ = iConfig.getUntrackedParameter<int>("zdc2DBins", 300);
  saveTrigger2D_ = iConfig.getUntrackedParameter<bool>("saveTrigger2D", true);
  zdcWideMax_  = iConfig.getUntrackedParameter<double>("zdcWideMax", 200000.);
}

void ZDCEnergyAnalyzer::beginJob()
{
  hNEvt = theOutputs->make<TH1D>("hNEvt", ";;events", 4, 0, 4);
  hNEvt->GetXaxis()->SetBinLabel(1, "processed");
  hNEvt->GetXaxis()->SetBinLabel(2, "centrality valid");
  hNEvt->GetXaxis()->SetBinLabel(3, "any configured trigger");
  hNEvt->GetXaxis()->SetBinLabel(4, "no configured trigger");

  const int nTrig = (int)triggerNames_.size();
  hTrigCount = theOutputs->make<TH1D>("hTrigCount", "fired configured paths;;events",
                                      nTrig > 0 ? nTrig : 1, 0, nTrig > 0 ? nTrig : 1);
  for(int i = 0; i < nTrig; i++)
    hTrigCount->GetXaxis()->SetBinLabel(i+1, triggerNames_[i].c_str());

  if(!dumpTriggerPattern_.empty())
  {
    // Alphanumeric axis: bins are created by name as paths are seen, so the
    // real menu comes out of the job without knowing it in advance.
    hAllTrig = theOutputs->make<TH1D>("hAllTrig", "every fired path matching the dump pattern;;events", 1, 0, 1);
    hAllTrig->SetCanExtend(TH1::kAllAxes);
  }

  hZdcPlus  = theOutputs->make<TH1D>("hZdcPlus",  ";ZDC+ energy sum;events", zdcHistBins_, zdcHistMin_, zdcHistMax_);
  hZdcMinus = theOutputs->make<TH1D>("hZdcMinus", ";ZDC- energy sum;events", zdcHistBins_, zdcHistMin_, zdcHistMax_);
  hZdcSum   = theOutputs->make<TH1D>("hZdcSum",   ";ZDC+ + ZDC- energy sum;events", zdcHistBins_, 2*zdcHistMin_, 2*zdcHistMax_);
  hZdcMinusVsPlus = theOutputs->make<TH2D>("hZdcMinusVsPlus", ";ZDC+ energy sum;ZDC- energy sum",
                                           zdc2DBins_, zdcHistMin_, zdcHistMax_, zdc2DBins_, zdcHistMin_, zdcHistMax_);

  hZdcPlusWide  = theOutputs->make<TH1D>("hZdcPlusWide",  "full range;ZDC+ energy sum;events", zdcWideBins_, zdcHistMin_, zdcWideMax_);
  hZdcMinusWide = theOutputs->make<TH1D>("hZdcMinusWide", "full range;ZDC- energy sum;events", zdcWideBins_, zdcHistMin_, zdcWideMax_);
  hZdcSumWide   = theOutputs->make<TH1D>("hZdcSumWide",   "full range;ZDC+ + ZDC- energy sum;events", zdcWideBins_, 2*zdcHistMin_, 2*zdcWideMax_);
  hZdcMinusVsPlusWide = theOutputs->make<TH2D>("hZdcMinusVsPlusWide", "full range;ZDC+ energy sum;ZDC- energy sum",
                                               zdc2DBins_, zdcHistMin_, zdcWideMax_, zdc2DBins_, zdcHistMin_, zdcWideMax_);
  hHF   = theOutputs->make<TH1D>("hHF",   ";#SigmaE^{HF}_{T} (GeV);events", 1240, 0, 6200);
  hNVtx = theOutputs->make<TH1D>("hNVtx", ";N_{vtx} (no cut applied);events", 51, -0.5, 50.5);

  hZdcPlusTrig.assign(nTrig, 0);
  hZdcMinusTrig.assign(nTrig, 0);
  hZdcSumTrig.assign(nTrig, 0);
  hZdcMinusVsPlusTrig.assign(nTrig, 0);
  for(int i = 0; i < nTrig; i++)
  {
    // Histogram names cannot carry the wildcard the path names end with.
    // Strip only a *trailing* "_v": a blanket ReplaceAll would also eat "_v"
    // occurring inside a path name (a hypothetical "..._veto..." would come
    // back as "...eto...").
    TString tag(triggerNames_[i].c_str());
    tag.ReplaceAll("*", "");
    if(tag.EndsWith("_v")) tag.Remove(tag.Length()-2);
    hZdcPlusTrig[i]  = theOutputs->make<TH1D>(Form("hZdcPlus_%s", tag.Data()),
        Form("%s;ZDC+ energy sum;events", triggerNames_[i].c_str()), zdcHistBins_, zdcHistMin_, zdcHistMax_);
    hZdcMinusTrig[i] = theOutputs->make<TH1D>(Form("hZdcMinus_%s", tag.Data()),
        Form("%s;ZDC- energy sum;events", triggerNames_[i].c_str()), zdcHistBins_, zdcHistMin_, zdcHistMax_);
    hZdcSumTrig[i]   = theOutputs->make<TH1D>(Form("hZdcSum_%s", tag.Data()),
        Form("%s;ZDC+ + ZDC- energy sum;events", triggerNames_[i].c_str()), zdcHistBins_, 2*zdcHistMin_, 2*zdcHistMax_);
    if(saveTrigger2D_)
      hZdcMinusVsPlusTrig[i] = theOutputs->make<TH2D>(Form("hZdcMinusVsPlus_%s", tag.Data()),
          Form("%s;ZDC+ energy sum;ZDC- energy sum", triggerNames_[i].c_str()),
          zdc2DBins_, zdcHistMin_, zdcHistMax_, zdc2DBins_, zdcHistMin_, zdcHistMax_);
  }

  if(saveTree_)
  {
    zdcTree = theOutputs->make<TTree>("zdcTree", "per-event ZDC sums and trigger decisions");
    zdcTree->Branch("run",      &runNb_,   "run/i");
    zdcTree->Branch("lumi",     &lsNb_,    "lumi/i");
    if(storeEvent_) zdcTree->Branch("event", &eventNb_, "event/l");
    zdcTree->Branch("zdcPlus",  &zdcPlus_, "zdcPlus/F");
    zdcTree->Branch("zdcMinus", &zdcMinus_,"zdcMinus/F");
    if(storeHF_) zdcTree->Branch("hfSumET", &hfSumET_, "hfSumET/F");
    if(doDigiSum_)
    {
      zdcTree->Branch("zdcPlusDigi",    &zdcPlusDigi_,    "zdcPlusDigi/F");
      zdcTree->Branch("zdcMinusDigi",   &zdcMinusDigi_,   "zdcMinusDigi/F");
      zdcTree->Branch("zdcPlusDigiId",  &zdcPlusDigiId_,  "zdcPlusDigiId/F");
      zdcTree->Branch("zdcMinusDigiId", &zdcMinusDigiId_, "zdcMinusDigiId/F");
      zdcTree->Branch("nZdcDigi",       &nZdcDigi_,       "nZdcDigi/I");
    }
    zdcTree->Branch("nPixel",   &nPixel_,  "nPixel/s");
    zdcTree->Branch("bx",       &bx_,      "bx/s");
    if(saveVertexInfo_)  zdcTree->Branch("nVtx",    &nVtx_,    "nVtx/I");
    if(useCentralityBin_) zdcTree->Branch("centBin", &centBin_, "centBin/I");
    // Bit i is triggerPathNames[i], in cfg order.
    zdcTree->Branch("trigBits", &trigBits_, "trigBits/l");
  }
}

int ZDCEnergyAnalyzer::findTriggerIndex(const edm::TriggerNames& names,
                                        const edm::TriggerResults& results,
                                        const std::string& path) const
{
  // Configured names carry a trailing wildcard ('..._v*'); strip it and match
  // the versioned menu entry by prefix.
  std::string stem = path;
  const size_t star = stem.find('*');
  if(star != std::string::npos) stem = stem.substr(0, star);

  std::vector<unsigned int> found;
  for(unsigned int i = 0; i < names.size(); i++)
    if(names.triggerName(i).find(stem) != std::string::npos && results.wasrun(i))
      found.push_back(i);

  if(found.empty()) return -1;
  // A menu can carry more than one version of a path; an accepted one wins.
  for(size_t k = 0; k < found.size(); k++)
    if(results.accept(found[k])) return (int)found[k];
  return (int)found[0];
}

void ZDCEnergyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  nProcessed_++;
  hNEvt->Fill(0.5);

  runNb_   = iEvent.id().run();
  bx_      = (unsigned short)iEvent.bunchCrossing();
  lsNb_    = iEvent.luminosityBlock();
  eventNb_ = iEvent.id().event();

  // ---- trigger decisions -------------------------------------------------
  trigBits_ = 0;
  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(token_triggerResults, triggerResults);
  if(triggerResults.isValid())
  {
    const edm::TriggerNames& names = iEvent.triggerNames(*triggerResults);
    for(size_t i = 0; i < triggerNames_.size(); i++)
    {
      const int idx = findTriggerIndex(names, *triggerResults, triggerNames_[i]);
      if(idx >= 0 && triggerResults->accept(idx))
      {
        trigBits_ |= (1ULL << i);
        hTrigCount->Fill((double)i + 0.5);
      }
    }
    if(hAllTrig)
      for(unsigned int i = 0; i < names.size(); i++)
        if(triggerResults->accept(i) &&
           names.triggerName(i).find(dumpTriggerPattern_) != std::string::npos)
          hAllTrig->Fill(names.triggerName(i).c_str(), 1.0);
  }
  hNEvt->Fill(trigBits_ ? 2.5 : 3.5);

  // ---- vertices: recorded, never required --------------------------------
  nVtx_ = -1;
  if(saveVertexInfo_)
  {
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(token_vertices, vertices);
    if(vertices.isValid())
    {
      nVtx_ = 0;
      for(reco::VertexCollection::const_iterator v = vertices->begin(); v != vertices->end(); ++v)
        if(!v->isFake() && (int)v->tracksSize() >= minVtxTracks_) nVtx_++;
    }
    hNVtx->Fill(nVtx_);
  }

  // ---- ZDC from digis (hardcoded calibration; the 2023 route) -------------
  zdcPlusDigi_ = -1.f; zdcMinusDigi_ = -1.f;
  zdcPlusDigiId_ = -1.f; zdcMinusDigiId_ = -1.f; nZdcDigi_ = 0;
  if(doDigiSum_)
  {
    edm::Handle<QIE10DigiCollection> zdcdigis;
    iEvent.getByToken(token_zdcdigi, zdcdigis);
    if(zdcdigis.isValid())
    {
      // (a) faithful port: channels picked by POSITION in the collection,
      //     EM minus 0-4, EM plus 12-16, HAD minus 8-11, HAD plus 20-23.
      // (b) the same sums selected by HcalZDCDetId instead.
      double pEM=0,pHD=0,mEM=0,mHD=0;          // positional
      double pEMi=0,pHDi=0,mEMi=0,mHDi=0;      // by DetId
      int idx=0;
      for(auto it = zdcdigis->begin(); it != zdcdigis->end(); ++it, ++idx)
      {
        const QIE10DataFrame digi = static_cast<const QIE10DataFrame>(*it);
        if(digi.samples() <= tsSig_ || tsPed_ < 0) continue;
        const double q =
            QWAna::ZDC2018::QIE10_regular_fC[ digi[tsSig_].adc() ][ digi[tsSig_].capid() ]
          - QWAna::ZDC2018::QIE10_regular_fC[ digi[tsPed_].adc() ][ digi[tsPed_].capid() ];
        nZdcDigi_++;

        if(idx >= 0  && idx < 5)  mEM += q;
        if(idx >= 8  && idx < 12) mHD += q;
        if(idx >= 12 && idx < 17) pEM += q;
        if(idx >= 20 && idx < 24) pHD += q;

        const HcalZDCDetId cell = digi.id();
        const int sec = cell.section(), z = cell.zside(), ch = cell.channel();
        if(sec == 1 && ch <= 5){ if(z > 0) pEMi += q; else mEMi += q; }
        if(sec == 2 && ch <= 4){ if(z > 0) pHDi += q; else mHDi += q; }
      }
      zdcPlusDigi_    = (Float_t)((pEM  * emWeight_ + pHD ) * calibPlus_);
      zdcMinusDigi_   = (Float_t)((mEM  * emWeight_ + mHD ) * calibMinus_);
      zdcPlusDigiId_  = (Float_t)((pEMi * emWeight_ + pHDi) * calibPlus_);
      zdcMinusDigiId_ = (Float_t)((mEMi * emWeight_ + mHDi) * calibMinus_);
      if(roundZDC_)
      {
        zdcPlusDigi_    = std::round(zdcPlusDigi_);
        zdcMinusDigi_   = std::round(zdcMinusDigi_);
        zdcPlusDigiId_  = std::round(zdcPlusDigiId_);
        zdcMinusDigiId_ = std::round(zdcMinusDigiId_);
      }
    }
  }

  // ---- ZDC ---------------------------------------------------------------
  zdcPlus_ = -1.f; zdcMinus_ = -1.f; hfSumET_ = -1.f; nPixel_ = 65535;
  edm::Handle<reco::Centrality> cent;
  iEvent.getByToken(token_centrality, cent);
  if(!cent.isValid())
  {
    // Without the centrality object there is no ZDC to record.  Count it and
    // keep going rather than aborting a long job.
    nNoCentrality_++;
    if(nNoCentrality_ <= 5)
      std::cout << "ZDCEnergyAnalyzer: centrality collection missing (run " << runNb_
                << ", ls " << lsNb_ << ")" << std::endl;
    // 2023 PromptReco has no centrality object at all.  With the digi sum on,
    // the event is still worth keeping -- that IS the measurement there.
    if(!(doDigiSum_ && !requireCentrality_)) return;
    if(zdcTree) zdcTree->Fill();
    return;
  }
  hNEvt->Fill(1.5);

  zdcPlus_  = cent->zdcSumPlus();
  zdcMinus_ = cent->zdcSumMinus();
  if(roundZDC_) { zdcPlus_ = std::round(zdcPlus_); zdcMinus_ = std::round(zdcMinus_); }
  hfSumET_  = cent->EtHFtowerSum();
  { const int np = (int)cent->multiplicityPixel();
    nPixel_ = np < 0 ? 0 : (np > 65534 ? 65534 : np); }

  centBin_ = -1;
  if(useCentralityBin_)
  {
    edm::Handle<int> cbin;
    iEvent.getByToken(token_centralitybin, cbin);
    if(cbin.isValid()) centBin_ = *cbin;
  }

  const double zdcSum = zdcPlus_ + zdcMinus_;
  hZdcPlus->Fill(zdcPlus_);
  hZdcMinus->Fill(zdcMinus_);
  hZdcSum->Fill(zdcSum);
  hZdcMinusVsPlus->Fill(zdcPlus_, zdcMinus_);
  hZdcPlusWide->Fill(zdcPlus_);
  hZdcMinusWide->Fill(zdcMinus_);
  hZdcSumWide->Fill(zdcSum);
  hZdcMinusVsPlusWide->Fill(zdcPlus_, zdcMinus_);
  hHF->Fill(hfSumET_);

  for(size_t i = 0; i < triggerNames_.size(); i++)
  {
    if(!(trigBits_ & (1ULL << i))) continue;
    hZdcPlusTrig[i]->Fill(zdcPlus_);
    hZdcMinusTrig[i]->Fill(zdcMinus_);
    hZdcSumTrig[i]->Fill(zdcSum);
    if(hZdcMinusVsPlusTrig[i]) hZdcMinusVsPlusTrig[i]->Fill(zdcPlus_, zdcMinus_);
  }

  if(zdcTree) zdcTree->Fill();
}

void ZDCEnergyAnalyzer::endJob()
{
  std::cout << "ZDCEnergyAnalyzer: " << nProcessed_ << " events processed";
  if(nNoCentrality_)
    std::cout << ", " << nNoCentrality_ << " without a centrality object (no ZDC recorded)";
  std::cout << std::endl;
}
