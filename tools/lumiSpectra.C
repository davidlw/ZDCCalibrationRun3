// ZDC spectra in bins of instantaneous luminosity, over a whole dataset.
//
// Saves the (L, ZDC) maps, which ARE the set of 1D spectra -- one projection per
// L bin -- so every later fit and parametrisation runs in seconds instead of
// re-reading a 12 GB tree.  One file per year.
//
// Binning choices that matter:
//   L     0.0001 wide below 0.003, 0.0002 above.  The curvature in mu(L) sits at
//         low L, and that is also where most of the integrated luminosity is.
//   ZDC   0-24000 in 50-unit bins: the n-neutron decomposition needs 1n, 2n and
//         3n, and in 2024 era A alone 2n already sits near 8000.
#include "TFile.h"
#include "TH2.h"
#include "TProfile.h"
#include "TString.h"
#include "ROOT/RDataFrame.hxx"
#include <set>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>

void lumiSpectra(const char* inFile, const char* csvFile, const char* outFile,
                 const char* year, int zbBit, int orBit,
                 const char* brP="zdcPlus", const char* brM="zdcMinus", const char* cut="1",
                 const char* jsonFile="")
{
  static std::map<ULong64_t,double> L; L.clear();
  { std::ifstream in(csvFile); std::string line;
    while(std::getline(in,line)){ if(line.empty()||line[0]=='#') continue;
      std::vector<std::string> p; std::stringstream ss(line); std::string t;
      while(std::getline(ss,t,',')) p.push_back(t); if(p.size()<7) continue;
      size_t c0=p[0].find(':'), c1=p[1].find(':'); if(c0==std::string::npos||c1==std::string::npos) continue;
      double d=0; try{ d=std::stod(p[5]); }catch(...){ continue; }
      L[(ULong64_t)std::stoul(p[0].substr(0,c0))*100000ULL+std::stoul(p[1].substr(0,c1))]=d/23.31; } }
  printf("%s: lumi table %zu LS\n",year,L.size()); fflush(stdout);

  // Optional certification JSON on top of the one already applied at production.
  // 2023 has a ZDC-specific one (Good_ZDC), which drops 31 of the 93 golden runs
  // -- including every run whose hardcoded digi sum reconstructs to nothing.
  static std::set<ULong64_t> good; good.clear();
  static bool useJson=false; useJson = (jsonFile && *jsonFile);
  if(useJson){
    std::ifstream in(jsonFile); std::stringstream ss; ss<<in.rdbuf(); std::string j=ss.str();
    size_t p=0;
    while((p=j.find('"',p))!=std::string::npos){
      size_t q=j.find('"',p+1); if(q==std::string::npos) break;
      const int run=atoi(j.substr(p+1,q-p-1).c_str());
      size_t a=j.find('[',q), b=j.find("]]",a); if(a==std::string::npos||b==std::string::npos) break;
      std::string body=j.substr(a+1,b-a);
      size_t r=0;
      while((r=body.find('[',r))!=std::string::npos){
        size_t c=body.find(',',r), e=body.find(']',c);
        const int lo=atoi(body.substr(r+1,c-r-1).c_str()), hi=atoi(body.substr(c+1,e-c-1).c_str());
        for(int l=lo;l<=hi;l++) good.insert((ULong64_t)run*100000ULL+l);
        r=e;
      }
      p=b+2;
    }
    printf("%s: certification JSON %zu lumi sections\n",year,good.size());
  }

  std::vector<double> edges;
  for(double x=0;      x<0.00299; x+=0.0001) edges.push_back(x);
  for(double x=0.0030; x<0.00801; x+=0.0002) edges.push_back(x);
  const int NL=(int)edges.size()-1;
  printf("%s: %d luminosity bins, %.4f to %.4f\n",year,NL,edges.front(),edges.back());

  ROOT::EnableImplicitMT(6);
  ROOT::RDataFrame d("zdc_ana/zdcTree", inFile);
  auto dd = d.Filter(cut)
             .Alias("P",brP).Alias("M",brM)
             .Define("Lv",[](UInt_t run, UInt_t ls){ auto it=L.find((ULong64_t)run*100000ULL+ls);
                                                     return it==L.end()? -1.0 : it->second; },{"run","lumi"})
             .Filter("Lv > 0")
             .Filter([](UInt_t run, UInt_t ls){ return !useJson || good.count((ULong64_t)run*100000ULL+ls)>0; },
                     {"run","lumi"});
  const int bits[2]={zbBit,orBit}; const char* nm[2]={"ZB","ZDCOR"};
  std::vector<ROOT::RDF::RResultPtr<TH2D>> hP,hM;
  std::vector<ROOT::RDF::RResultPtr<TProfile>> pL;
  std::vector<ROOT::RDF::RResultPtr<ULong64_t>> cnt;
  for(int i=0;i<2;i++){
    const int b=bits[i];
    auto f=dd.Filter([b](ULong64_t tb){ return (bool)((tb>>b)&1ULL); },{"trigBits"});
    hP.push_back(f.Histo2D({Form("hP_%s_%s",year,nm[i]),Form("%s %s;inst. lumi (/#mub/s);ZDC+ (a.u.)",year,nm[i]),
                            NL,&edges[0],480,0,24000},"Lv","P"));
    hM.push_back(f.Histo2D({Form("hM_%s_%s",year,nm[i]),Form("%s %s;inst. lumi (/#mub/s);ZDC- (a.u.)",year,nm[i]),
                            NL,&edges[0],480,0,24000},"Lv","M"));
    // the mean L inside each bin, which is not the bin centre when the
    // luminosity profile is steep
    pL.push_back(f.Profile1D({Form("pL_%s_%s",year,nm[i]),Form("%s %s;inst. lumi bin;<L>",year,nm[i]),
                              NL,&edges[0]},"Lv","Lv"));
    cnt.push_back(f.Count());
  }
  auto nAll=dd.Count();
  TFile* fo=new TFile(outFile,"RECREATE","",505);
  for(int i=0;i<2;i++){ hP[i]->Write(); hM[i]->Write(); pL[i]->Write(); }
  fo->Close();
  printf("%s: %llu events with a luminosity match\n",year,*nAll);
  for(int i=0;i<2;i++) printf("   bit %d %-6s %12llu\n",bits[i],nm[i],*cnt[i]);
  fflush(stdout);
}
