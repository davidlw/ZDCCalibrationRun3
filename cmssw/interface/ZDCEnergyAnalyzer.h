#ifndef ZDCEnergyAnalyzer_
#define ZDCEnergyAnalyzer_

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDetId/interface/HcalZDCDetId.h"

class TH1D;
class TH2D;
class TTree;

// Total ZDC energy per side in zero-bias and ZDC-OR triggered events.
//
// Deliberately does NOT require a reconstructed vertex.  SoSMeanPtAnalyzer
// returns early on nVertices<1, where nVertices counts only non-fake vertices
// with >=2 tracks, so it never sees the events this study is about: a 0n0n UPC
// event, or a zero-bias event with no tracks, has no such vertex.  Nothing here
// is measured with respect to a vertex, so there is nothing to require.
//
// ZDC comes straight off reco::Centrality (zdcSumPlus/zdcSumMinus), which is
// the calibrated sum for 2024-2026 data.  No digi unpacking, no electronics
// map, no QWzdcreco -- that machinery is only needed for 2023.
class ZDCEnergyAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {

 public:
  explicit ZDCEnergyAnalyzer(const edm::ParameterSet&);

 private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // Substring match against the versioned HLT names, the same way
  // EventInfoTreeProducer resolves its configured paths.  Returns -1 if the
  // path is absent from the menu.
  int findTriggerIndex(const edm::TriggerNames&, const edm::TriggerResults&,
                       const std::string&) const;

  edm::Service<TFileService> theOutputs;

  edm::EDGetTokenT<reco::Centrality> token_centrality;
  edm::EDGetTokenT<int> token_centralitybin;
  edm::EDGetTokenT<edm::TriggerResults> token_triggerResults;
  edm::EDGetTokenT<reco::VertexCollection> token_vertices;
  edm::EDGetTokenT<QIE10DigiCollection> token_zdcdigi;

  // Configured HLT paths, in the order their bits are packed into trigBits_.
  std::vector<std::string> triggerNames_;
  bool useCentralityBin_;
  bool saveTree_;
  bool saveVertexInfo_;
  int  minVtxTracks_;
  // Non-empty: label-fill every fired path whose name contains this string into
  // hAllTrig, so the menu can be read off a short test run.  See the cfg.
  std::string dumpTriggerPattern_;

  // Narrow range: the 0n/1n region, where essentially all the events are.
  // Starts below zero on purpose so pedestal-negative values are visible
  // instead of being buried in the underflow bin.
  int    zdcHistBins_;
  double zdcHistMin_;
  double zdcHistMax_;
  // Wide range: the Xn tail, inclusive only -- per-trigger tails are better
  // taken from the tree.
  int    zdcWideBins_;
  double zdcWideMax_;
  int    zdc2DBins_;
  bool   saveTrigger2D_;

  // ---- hardcoded digi-level ZDC sum (the 2023 route) ----------------------
  // 2023 has no usable ZDC in reco::Centrality: PromptReco has no centrality
  // object at all, and the 14Feb2025 rereco has one whose sums show no neutron
  // peak structure.  So 2023 must be summed straight from the QIE10 digis with
  // a hardcoded calibration, as in the 2023 SoSMeanPtAnalyzer.
  //
  // Signal is a plain TS(sig) - TS(ped) difference: no conditions pedestal and
  // no template fit, so there is no out-of-time subtraction at all.  Note this
  // uses a LOCAL baseline, which should largely cancel the rate-proportional
  // pedestal that drives the additive term in the Run 3 route (see the study
  // summary, sections 16-17) -- a testable difference between the two methods.
  bool   doDigiSum_;
  bool   requireCentrality_;   // off for 2023, where centrality may be absent
  double emWeight_;            // EM weight relative to HAD (0.1)
  double calibPlus_, calibMinus_;
  int    tsSig_, tsPed_;
  bool   storeErrors_;

  // ---- tree slimming -------------------------------------------------------
  // Measured on the 2025 tree: 15.5 B/event, of which the event number (21%)
  // and hfSumET (22%) were never read by any macro, and the two ZDC floats
  // (42%) carry a random mantissa that defeats compression.  Rounding the sums
  // to integers, dropping those two branches, storing nPixel/bx as 16-bit and
  // writing with ZSTD gives 4.8 B/event, i.e. 3.2x smaller, with every reader
  // unchanged (the ZDC branches stay Float_t).
  bool   storeEvent_;          // event number branch (off)
  bool   storeHF_;             // hfSumET branch (off)
  bool   roundZDC_;            // round ZDC sums to integers (on)
  int    compression_;         // ROOT compression setting, alg*100+level

  // monitoring histograms -- the tree is what the analysis is done from
  TH1D* hNEvt;
  TH1D* hTrigCount;
  TH1D* hAllTrig;
  TH1D* hZdcPlus;
  TH1D* hZdcMinus;
  TH1D* hZdcSum;
  TH2D* hZdcMinusVsPlus;
  TH1D* hZdcPlusWide;
  TH1D* hZdcMinusWide;
  TH1D* hZdcSumWide;
  TH2D* hZdcMinusVsPlusWide;
  TH1D* hHF;
  TH1D* hNVtx;
  // per configured trigger
  std::vector<TH1D*> hZdcPlusTrig;
  std::vector<TH1D*> hZdcMinusTrig;
  std::vector<TH1D*> hZdcSumTrig;
  std::vector<TH2D*> hZdcMinusVsPlusTrig;

  TTree* zdcTree;
  // Two independent channel selections, stored side by side so the fragile one
  // can be validated: *Digi uses POSITION in the digi collection exactly as the
  // original 2023 code does, *DigiId selects by HcalZDCDetId instead.  They
  // should agree; if they diverge the collection ordering assumption is wrong.
  Float_t zdcPlusDigi_,   zdcMinusDigi_;
  Float_t zdcPlusDigiId_, zdcMinusDigiId_;
  Int_t   nZdcDigi_;
  unsigned int   runNb_;
  unsigned int   lsNb_;
  unsigned long long eventNb_;
  float          zdcPlus_;
  float          zdcMinus_;
  float          hfSumET_;
  // 16-bit: nPixel saturates at 65535 (the UPC paths cap it at 10000 anyway),
  // bx is at most 3564.  65535 in nPixel marks "no centrality object".
  unsigned short nPixel_;
  // Bunch crossing number.  The ZDC integrates over a window that catches
  // neighbouring crossings, so out-of-time pileup depends on position within
  // the bunch train, which no in-time proxy (nPixel, nVtx) can reproduce.
  unsigned short bx_;
  int            nVtx_;
  int            centBin_;
  unsigned long long trigBits_;

  unsigned long long nProcessed_;
  unsigned long long nNoCentrality_;
};

#endif  // ZDCEnergyAnalyzer_
