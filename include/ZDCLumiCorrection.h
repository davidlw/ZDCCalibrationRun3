#ifndef ZDCLumiCorrection_h
#define ZDCLumiCorrection_h
//
// Luminosity correction for the ZDC energy sum, Run 3 PbPb (2023-2026).
//
//   ZDC_corr = ZDC_raw * muRef / mu1(L)
//
// The ZDC energy scale rises with instantaneous luminosity: the 1n peak moves by
// up to +57% across a fill.  mu1(L) is the measured 1n peak position, fitted
// with a cubic per year and per side; dividing by it and rescaling to a common
// reference puts the scale back where it belongs at every luminosity.
//
// This is a SINGLE-PARAMETER correction, tuned on the 1n peak.  It is what an
// Xn (>=1 neutron) selection needs: the peak lands at a fixed place and the
// valley below it moves with the peak, so a cut at a fixed fraction of muRef
// stays in the valley at every luminosity.  It does NOT make the peaks correctly
// SPACED -- for that you need the two-component form with b(L) and k(L), which
// matters only if you are counting neutrons (1n vs 2n vs 3n).
//
// The per-event luminosity is not in MINIAOD (LumiScalers is empty, there is no
// OnlineLuminosityRecord), so it comes from a brilcalc --byls table joined on
// (run, lumiSection).  Make that table with
//
//   brilcalc lumi -c web -i <golden json> --byls -u /ub -o lumi_byls_<YEAR>HI.csv
//
// and note that field 5 is the luminosity DELIVERED PER LUMI SECTION; one lumi
// section is 23.31 s, so L = field5 / 23.31 in /ub/s.  That division is done
// here, so pass the csv exactly as brilcalc writes it.
//
// Usage, offline:
//     ZDCLumiCorrection zdcCorr;
//     zdcCorr.addLumiTable("lumi_byls_2025HI.csv");
//     ...
//     double zp = zdcCorr.correct(run, lumi, ZDCLumiCorrection::Plus,  zdcPlusRaw);
//     if (!zdcCorr.valid()) { /* no luminosity for this lumi section */ }
//
// Usage, in an EDProducer: call addLumiTable() once in the constructor, and
// correct() per event.  Looking L up per lumi block rather than per event is
// cheaper, but the map lookup here is already negligible.
//
// Honest limits, measured:
//   * The cubic reproduces the measured mu1(L) to about 1% (2024, 2026), 2%
//     (2025) and 1-2% (2023, where that is consistent with statistics alone).
//   * Outside the fitted luminosity range L is CLAMPED to the range edge rather
//     than extrapolated: a cubic extrapolates violently.
//   * 2023 barely needs this at all -- its scale moves ~1-2% in total, because
//     its TS2-TS1 estimator uses a local baseline.
//   * 2024 before run 388556 is NOT supported: that period contains three
//     calibration states with ~8% run-to-run scatter at fixed luminosity, and
//     needs a per-run normalisation first.
//
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

class ZDCLumiCorrection {
public:
  enum Side { Plus = 0, Minus = 1 };

  ZDCLumiCorrection() : ok_(false) {}

  // Load a brilcalc --byls csv.  Call once per year you intend to correct; the
  // tables are keyed by run so several can be loaded together.
  bool addLumiTable(const std::string& csvFile) {
    std::ifstream in(csvFile.c_str());
    if (!in) { std::fprintf(stderr, "ZDCLumiCorrection: cannot open %s\n", csvFile.c_str()); return false; }
    std::string line; size_t n = 0;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::vector<std::string> f; std::stringstream ss(line); std::string t;
      while (std::getline(ss, t, ',')) f.push_back(t);
      if (f.size() < 7) continue;
      const size_t c0 = f[0].find(':'), c1 = f[1].find(':');
      if (c0 == std::string::npos || c1 == std::string::npos) continue;
      double delivered = 0;
      try { delivered = std::stod(f[5]); } catch (...) { continue; }
      const unsigned run = std::strtoul(f[0].substr(0, c0).c_str(), 0, 10);
      const unsigned ls  = std::strtoul(f[1].substr(0, c1).c_str(), 0, 10);
      lumi_[key(run, ls)] = delivered / 23.31;   // delivered per LS -> /ub/s
      ++n;
    }
    std::printf("ZDCLumiCorrection: %s -> %zu lumi sections\n", csvFile.c_str(), n);
    return n > 0;
  }

  // Instantaneous luminosity for this lumi section, or -1 if it is not in the
  // table (a non-certified section, or a table for the wrong year).
  double instLumi(unsigned run, unsigned lumiSection) const {
    std::map<unsigned long long, double>::const_iterator it = lumi_.find(key(run, lumiSection));
    return it == lumi_.end() ? -1.0 : it->second;
  }

  // The fitted 1n peak position at this luminosity, in raw ZDC units.
  double mu1(unsigned run, Side side, double L) const {
    const Era* e = era(run);
    if (!e) return -1.0;
    const double x = std::min(std::max(L, e->lmin), e->lmax);   // clamp, never extrapolate
    const double* p = e->p[side];
    return p[0] + x * (p[1] + x * (p[2] + x * p[3]));
  }

  // The correction factor alone, if you would rather apply it yourself.
  double factor(unsigned run, unsigned lumiSection, Side side) const {
    ok_ = false;
    const Era* e = era(run);
    if (!e) return 1.0;
    const double L = instLumi(run, lumiSection);
    if (L <= 0) return 1.0;
    const double m = mu1(run, side, L);
    if (m <= 0) return 1.0;
    ok_ = true;
    return e->muRef[side] / m;
  }

  // Corrected ZDC sum.  Check valid() afterwards: when the lumi section has no
  // luminosity, or the run is outside the supported ranges, the raw value is
  // returned unchanged rather than silently scaled by something wrong.
  double correct(unsigned run, unsigned lumiSection, Side side, double zdcRaw) const {
    return zdcRaw * factor(run, lumiSection, side);
  }

  // True if the last correct()/factor() call actually applied a correction.
  bool valid() const { return ok_; }

  // The reference scale a corrected value is expressed in (the fitted mu1 at
  // L = 0.003 /ub/s).  An Xn cut should be quoted as a fraction of this.
  double muRef(unsigned run, Side side) const {
    const Era* e = era(run);
    return e ? e->muRef[side] : -1.0;
  }

private:
  struct Era {
    const char* name;
    unsigned runMin, runMax;
    double lmin, lmax;        // fitted luminosity range; L is clamped to it
    double p[2][4];           // cubic coefficients, [side][order]
    double muRef[2];          // mu1 at L = 0.003 /ub/s
  };

  static unsigned long long key(unsigned run, unsigned ls) {
    return (unsigned long long)run * 100000ULL + ls;
  }

  const Era* era(unsigned run) const {
    for (size_t i = 0; i < nEras(); ++i)
      if (run >= eras()[i].runMin && run <= eras()[i].runMax) return &eras()[i];
    return 0;
  }

  // Fitted on the full primary datasets, ZDC1nOR trigger, golden JSON (2023 also
  // with the ZDC-specific certification, Cert_..._Good_ZDC_Golden.json).
  static const Era* eras() {
    static const Era e[] = {
      // 2023.  Scale is nearly flat in L; included for completeness.
      { "2023",      374288, 375823, 0.0005, 0.0035,
        { { 2.549355e+03, -3.636667e+05,  1.841071e+08, -2.631790e+10 },
          { 2.910405e+03, -6.746395e+05,  3.242357e+08, -4.458352e+10 } },
        { 2404.74, 2600.85 } },
      // 2024, runs >= 388556 only.  Earlier 2024 has three calibration states
      // and needs a per-run normalisation before this can be used.
      { "2024final", 388556, 388784, 0.0010, 0.0060,
        { { 2.170344e+03,  7.848430e+04,  1.705695e+07, -1.261845e+09 },
          { 2.097376e+03,  1.112956e+05,  1.877651e+07, -4.344867e+08 } },
        { 2525.24, 2588.52 } },
      { "2025",      399465, 400426, 0.0010, 0.0060,
        { { 2.242950e+03,  6.227552e+04,  4.227594e+07, -4.306951e+09 },
          { 2.095270e+03,  5.765722e+04,  5.397045e+07, -4.317196e+09 } },
        { 2693.97, 2637.41 } },
      // 2026.  ZDC+ saturates above L ~ 0.0045; the cubic follows it, but treat
      // ZDC+ above that luminosity with care.
      { "2026",      404469, 404926, 0.0010, 0.0060,
        { { 2.095278e+03,  1.086939e+05,  2.665948e+07, -4.102765e+09 },
          { 1.923198e+03,  9.078184e+04,  3.040025e+07, -1.936089e+09 } },
        { 2550.52, 2416.87 } },
    };
    return e;
  }
  static size_t nEras() { return 4; }

  std::map<unsigned long long, double> lumi_;
  mutable bool ok_;
};

#endif
