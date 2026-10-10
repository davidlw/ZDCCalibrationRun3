// Run-level closure: derive the correction on ONE long run, apply it to a
// DIFFERENT long run, and ask whether the second run's 1n peak becomes
// luminosity independent.  Deriving and applying on the same data proves
// nothing; two runs from different fills is the honest version.
//
// Per luminosity bin the 1n and 2n peaks give the two components directly:
//     mu_n = b + n*k   ->   k = mu_2 - mu_1,   b = 2*mu_1 - mu_2
// b(L) and k(L) are then fitted linear in L on the derivation run, and the
// correction (ZDC - b(L)) * k_ref/k(L) with k_ref = k(L->0) is applied to the
// other run's measured peak.
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TString.h"
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>

static const char* FIGD2 = "/Users/wl33/Library/CloudStorage/Dropbox/UPC/ZDC/figures";

static std::map<ULong64_t,double> loadL(const char* csv){
  std::map<ULong64_t,double> L; std::ifstream in(csv); std::string line;
  while(std::getline(in,line)){ if(line.empty()||line[0]=='#') continue;
    std::vector<std::string> p; std::stringstream ss(line); std::string t;
    while(std::getline(ss,t,',')) p.push_back(t); if(p.size()<7) continue;
    size_t c0=p[0].find(':'), c1=p[1].find(':'); if(c0==std::string::npos||c1==std::string::npos) continue;
    double d=0; try{ d=std::stod(p[5]); }catch(...){ continue; }
    L[(ULong64_t)std::stoul(p[0].substr(0,c0))*100000ULL+std::stoul(p[1].substr(0,c1))]=d/23.31; }
  return L;
}

static TH2D* mapLvsZ(const char* file, int run, ULong64_t mask, const char* br,
                     std::map<ULong64_t,double>& L, const char* nm){
  TFile* f=TFile::Open(file); if(!f||f->IsZombie()) return 0;
  TTree* t=(TTree*)f->Get("zdc_ana/zdcTree"); if(!t) return 0;
  TH2D* h=new TH2D(nm,"",34,0,0.0068,320,0,24000);
  h->SetDirectory(0);
  TTreeReader r(t);
  TTreeReaderValue<UInt_t> rn(r,"run"), ls(r,"lumi");
  TTreeReaderValue<ULong64_t> tb(r,"trigBits");
  TTreeReaderValue<Float_t> z(r,br);
  while(r.Next()){
    if((int)*rn!=run || !(*tb & mask)) continue;
    auto it=L.find((ULong64_t)(*rn)*100000ULL+*ls);
    if(it==L.end()||it->second<=0) continue;
    h->Fill(it->second,*z);
  }
  return h;
}

static double peakNear(TH1* h, double seed, double win){
  const double w=win*seed;
  int b1=h->FindBin(seed-w), b2=h->FindBin(seed+w), bm=0; double mx=0;
  for(int b=std::max(b1,1);b<=std::min(b2,h->GetNbinsX());b++)
    if(h->GetBinContent(b)>mx){ mx=h->GetBinContent(b); bm=b; }
  return bm? h->GetBinCenter(bm) : 0;
}

// 1n and 2n peak per luminosity bin
static int peaks12(TH2D* h, std::vector<double>& x, std::vector<double>& m1,
                   std::vector<double>& m2, double nmin=4000){
  x.clear(); m1.clear(); m2.clear();
  double seed=0;
  { TH1D* a=h->ProjectionY("a_all");
    int bm=0; double v=0;
    for(int b=a->FindBin(1400.);b<=a->GetNbinsX();b++) if(a->GetBinContent(b)>v){v=a->GetBinContent(b);bm=b;}
    seed=bm? a->GetBinCenter(bm):2500; }
  for(int bx=1;bx<=h->GetNbinsX();bx++){
    TH1D* c=h->ProjectionY(Form("%s_c%d",h->GetName(),bx),bx,bx);
    if(c->Integral(c->FindBin(1200.),c->GetNbinsX())<nmin) continue;
    const double p1=peakNear(c,seed,0.22); if(p1<=0) continue;
    TF1 g1("g1","gaus",p1*0.75,p1*1.28);
    g1.SetParameters(c->GetBinContent(c->FindBin(p1)),p1,0.2*p1);
    g1.SetParLimits(1,p1*0.82,p1*1.18); g1.SetParLimits(2,0.05*p1,0.40*p1);
    if(c->Fit(&g1,"QNR")!=0) continue;
    const double mu1=g1.GetParameter(1);
    // 2n sits near 2*mu1; search a narrow window so a shoulder cannot pull it
    const double p2=peakNear(c,2*mu1,0.16); if(p2<=0) continue;
    TF1 g2("g2","gaus",p2*0.86,p2*1.14);
    g2.SetParameters(c->GetBinContent(c->FindBin(p2)),p2,0.15*p2);
    g2.SetParLimits(1,2*mu1*0.86,2*mu1*1.14); g2.SetParLimits(2,0.05*p2,0.30*p2);
    if(c->Fit(&g2,"QNR")!=0) continue;
    const double mu2=g2.GetParameter(1);
    if(mu2<1.4*mu1 || mu2>2.6*mu1) continue;          // not a credible 2n
    seed=mu1;
    x.push_back(h->GetXaxis()->GetBinCenter(bx)); m1.push_back(mu1); m2.push_back(mu2);
  }
  return (int)x.size();
}

void closure(const char* year, const char* fileDer, int runDer,
             const char* fileApp, int runApp, const char* csv,
             ULong64_t mask, const char* brP, const char* brM)
{
  gStyle->SetOptStat(0); gStyle->SetTitleFontSize(0.048);
  gStyle->SetTitleAlign(23); gStyle->SetTitleX(0.56); gStyle->SetTitleY(0.985);
  gSystem->mkdir(Form("%s/closure",FIGD2),kTRUE);
  auto L=loadL(csv);
  TCanvas* c=new TCanvas(Form("ccl_%s",year),"",1500,660); c->Divide(2,1,0.003,0.003);
  for(int s=0;s<2;s++){
    const char* br = s? brM : brP;
    TH2D* hD=mapLvsZ(fileDer,runDer,mask,br,L,Form("hD_%s_%d",year,s));
    TH2D* hA=mapLvsZ(fileApp,runApp,mask,br,L,Form("hA_%s_%d",year,s));
    if(!hD||!hA) { printf("%s: missing maps\n",year); return; }
    std::vector<double> xD,m1D,m2D, xA,m1A,m2A;
    const int nD=peaks12(hD,xD,m1D,m2D), nA=peaks12(hA,xA,m1A,m2A);
    if(nD<4||nA<3){ printf("%s ZDC%s: too few bins (derive %d, apply %d)\n",year,s?"-":"+",nD,nA); continue; }
    // b(L), k(L) from the derivation run
    TGraph gb,gk;
    for(int i=0;i<nD;i++){ gb.SetPoint(i,xD[i],2*m1D[i]-m2D[i]); gk.SetPoint(i,xD[i],m2D[i]-m1D[i]); }
    // b(L) and k(L) linear forces mu_1n linear, but mu_1n is measurably curved:
    // ZDC- is convex (slope 33-115% steeper at high L) in every year, so the
    // end of the fill moves far less per unit L than a straight line predicts.
    // Fit both components quadratic and compare against the linear form.
    TF1 fb("fb","pol1",0,0.0068), fk("fk","pol1",0,0.0068);
    gb.Fit(&fb,"QN"); gk.Fit(&fk,"QN");
    TF1 fb2("fb2","pol2",0,0.0068), fk2("fk2","pol2",0,0.0068);
    gb.Fit(&fb2,"QN"); gk.Fit(&fk2,"QN");
    const double kref2=fk2.Eval(0.);
    const double kref=fk.Eval(0.);
    // apply to the OTHER run's measured 1n peak
    TGraph *graw=new TGraph(), *gcor=new TGraph(), *gcor2=new TGraph();
    for(int i=0;i<nA;i++){
      const double Lv=xA[i];
      graw->SetPoint(i,Lv,m1A[i]);
      gcor->SetPoint(i,Lv,(m1A[i]-fb.Eval(Lv))*kref/fk.Eval(Lv));
      gcor2->SetPoint(i,Lv,(m1A[i]-fb2.Eval(Lv))*kref2/fk2.Eval(Lv));
    }
    // spread of the corrected points, a fairer measure than first-vs-last
    auto spread=[&](TGraph* g){ double lo=1e9,hi=-1e9,sum=0;
      for(int i=0;i<nA;i++){ lo=std::min(lo,g->GetY()[i]); hi=std::max(hi,g->GetY()[i]); sum+=g->GetY()[i]; }
      return 100.*(hi-lo)/(sum/nA); };
    const double r0=graw->GetY()[0], r1=graw->GetY()[nA-1];
    const double c0=gcor->GetY()[0], c1=gcor->GetY()[nA-1];
    TPad* p=(TPad*)c->cd(s+1); p->SetLeftMargin(0.14); p->SetRightMargin(0.04); p->SetTopMargin(0.14);
    graw->SetMarkerStyle(20); graw->SetMarkerColor(kGray+2); graw->SetLineColor(kGray+2); graw->SetLineWidth(2);
    gcor->SetMarkerStyle(21); gcor->SetMarkerColor(kAzure+2); gcor->SetLineColor(kAzure+2); gcor->SetLineWidth(2);
    double lo=1e9,hi=-1e9;
    for(int i=0;i<nA;i++){ lo=std::min({lo,graw->GetY()[i],gcor->GetY()[i]}); hi=std::max({hi,graw->GetY()[i],gcor->GetY()[i]}); }
    graw->SetTitle(Form("#splitline{%s   ZDC%s   derived on run %d, applied to run %d}{}",year,s?"-":"+",runDer,runApp));
    graw->GetXaxis()->SetTitle("inst. lumi (/#mub/s)");
    graw->GetYaxis()->SetTitle("1n peak (a.u.)"); graw->GetYaxis()->SetTitleOffset(1.5);
    graw->GetXaxis()->SetLimits(0,0.0070); graw->SetMinimum(lo*0.88); graw->SetMaximum(hi*1.15);
    graw->Draw("APL"); gcor->Draw("PL same");
    TLegend* lg=new TLegend(0.42,0.70,0.96,0.845); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.033);
    lg->AddEntry(graw,Form("raw        %+.1f%% over the range",100*(r1-r0)/r0),"pl");
    lg->AddEntry(gcor,Form("corrected  %+.1f%%",100*(c1-c0)/c0),"pl");
    lg->Draw();
    const double q0=gcor2->GetY()[0], q1=gcor2->GetY()[nA-1];
    printf("%-6s ZDC%s  run %d -> %d:  raw %+.1f%%   corrected(linear) %+.1f%% [spread %.1f%%]   corrected(quadratic) %+.1f%% [spread %.1f%%]\n",
           year,s?"-":"+",runDer,runApp,100*(r1-r0)/r0,100*(c1-c0)/c0,spread(gcor),100*(q1-q0)/q0,spread(gcor2));
  }
  c->SaveAs(Form("%s/closure/closure_%s.png",FIGD2,year));
  delete c;
}
