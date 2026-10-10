// b(L) and k(L) from the luminosity-binned spectra.
//
// Per L bin the 1n, 2n and 3n peaks give the two components directly:
//     mu_n = b + n*k
// Fitting 1n and 2n gives k = mu_2 - mu_1 and b = 2*mu_1 - mu_2; the 3n peak,
// where it is resolved, is then a consistency check on the linearity in n and
// is NOT used to define b or k.
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include "TString.h"
#include <vector>
#include <algorithm>

static const char* FIGD3="/Users/wl33/Library/CloudStorage/Dropbox/UPC/ZDC/figures";

static double peakNear(TH1* h, double seed, double win){
  const double w=win*seed;
  int b1=h->FindBin(seed-w), b2=h->FindBin(seed+w), bm=0; double mx=0;
  for(int b=std::max(b1,1);b<=std::min(b2,h->GetNbinsX());b++)
    if(h->GetBinContent(b)>mx){ mx=h->GetBinContent(b); bm=b; }
  return bm? h->GetBinCenter(bm) : 0;
}
static bool fitPeak(TH1* h, double seed, double win, double lo, double hi,
                    double& mu, double& sg){
  const double p=peakNear(h,seed,win); if(p<=0) return false;
  TF1 g("gp","gaus",p*lo,p*hi);
  g.SetParameters(h->GetBinContent(h->FindBin(p)),p,0.18*p);
  g.SetParLimits(1,p*(1-0.8*win),p*(1+0.8*win));
  g.SetParLimits(2,0.04*p,0.42*p);
  if(h->Fit(&g,"QNR")!=0) return false;
  mu=g.GetParameter(1); sg=fabs(g.GetParameter(2));
  return mu>500 && mu<40000 && sg/mu<0.5;
}

void fitLumiDep(const char* year, const char* trig="ZDCOR", double nmin=20000){
  gStyle->SetOptStat(0); gStyle->SetTitleFontSize(0.050);
  gStyle->SetTitleAlign(23); gStyle->SetTitleX(0.56); gStyle->SetTitleY(0.995);
  gSystem->mkdir(Form("%s/lumi_dep",FIGD3),kTRUE);
  TFile* f=TFile::Open(Form("lumispectra_%s.root",year));
  if(!f||f->IsZombie()){ printf("no spectra for %s\n",year); return; }
  TProfile* pL=(TProfile*)f->Get(Form("pL_%s_%s",year,trig));
  const char* sn[2]={"+","-"};
  std::vector<double> x[2],mu1[2],mu2[2],mu3[2],kk[2],bb[2],res[2];
  for(int s=0;s<2;s++){
    TH2D* h=(TH2D*)f->Get(Form("h%c_%s_%s",s?'M':'P',year,trig));
    if(!h){ printf("missing map %s %s\n",year,trig); return; }
    double seed=0;
    { TH1D* a=h->ProjectionY("aall");
      int bm=0; double v=0;
      for(int b=a->FindBin(1400.);b<=a->GetNbinsX();b++) if(a->GetBinContent(b)>v){v=a->GetBinContent(b);bm=b;}
      seed=bm? a->GetBinCenter(bm):2500; }
    for(int bx=1;bx<=h->GetNbinsX();bx++){
      TH1D* c=h->ProjectionY(Form("p%s%d%d",year,s,bx),bx,bx);
      if(c->Integral(c->FindBin(1200.),c->GetNbinsX())<nmin) continue;
      double m1,s1,m2,s2,m3,s3;
      if(!fitPeak(c,seed,0.22,0.74,1.30,m1,s1)) continue;
      if(!fitPeak(c,2*m1,0.16,0.87,1.13,m2,s2)) continue;
      if(m2<1.5*m1 || m2>2.5*m1) continue;
      const bool ok3 = fitPeak(c,3*m1,0.12,0.90,1.10,m3,s3) && m3>2.4*m1 && m3<3.6*m1;
      seed=m1;
      const double Lv = pL->GetBinEntries(bx)>0? pL->GetBinContent(bx) : h->GetXaxis()->GetBinCenter(bx);
      x[s].push_back(Lv); mu1[s].push_back(m1); mu2[s].push_back(m2);
      mu3[s].push_back(ok3? m3 : 0);
      kk[s].push_back(m2-m1); bb[s].push_back(2*m1-m2); res[s].push_back(s1/m1);
    }
    printf("\n%s %s ZDC%s : %zu luminosity bins\n",year,trig,sn[s],x[s].size());
    printf("   %8s %8s %8s %8s %8s %8s %7s\n","L","mu1","mu2","mu3","k","b","sig/mu");
    for(size_t i=0;i<x[s].size();i+=3)
      printf("   %8.4f %8.0f %8.0f %8.0f %8.0f %8.0f %7.3f\n",
             x[s][i],mu1[s][i],mu2[s][i],mu3[s][i],kk[s][i],bb[s][i],res[s][i]);
    // 3n consistency: mu3 should equal b + 3k where it is resolved
    int n3=0; double d3=0;
    for(size_t i=0;i<x[s].size();i++) if(mu3[s][i]>0){ n3++; d3+=fabs(mu3[s][i]-(bb[s][i]+3*kk[s][i]))/mu3[s][i]; }
    if(n3) printf("   3n resolved in %d bins; |mu3 - (b+3k)|/mu3 = %.1f%% on average\n",n3,100*d3/n3);
  }
  // Two figures: the measured quantities (the 1n peak and its resolution) and
  // the model parameters derived from the peak spacing (k and b).
  const char* ttl[4]={"1n peak #mu_{1} (a.u.)","per-neutron scale k (a.u.)","offset b (a.u.)","resolution #sigma/#mu"};
  const int col[2]={kAzure+2,kRed+1};
  const int pairs[2][2]={{0,3},{1,2}};                    // {mu1, sigma/mu} and {k, b}
  const char* nick[2]={"peak","bk"};
  for(int fig=0;fig<2;fig++){
    TCanvas* c=new TCanvas(Form("cld_%s_%d",year,fig),"",1500,660); c->Divide(2,1,0.004,0.004);
    for(int j=0;j<2;j++){
      const int q=pairs[fig][j];
      TPad* p=(TPad*)c->cd(j+1); p->SetLeftMargin(0.13); p->SetRightMargin(0.04); p->SetTopMargin(0.085);
      TLegend* lg=new TLegend(0.62,0.74,0.95,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.040);
      double lo=1e9,hi=-1e9;
      for(int s=0;s<2;s++) for(size_t i=0;i<x[s].size();i++){
        const double v = q==0? mu1[s][i] : q==1? kk[s][i] : q==2? bb[s][i] : res[s][i];
        lo=std::min(lo,v); hi=std::max(hi,v); }
      if(lo>hi) continue;
      for(int s=0;s<2;s++){
        TGraph* g=new TGraph();
        for(size_t i=0;i<x[s].size();i++)
          g->SetPoint(g->GetN(),x[s][i], q==0? mu1[s][i] : q==1? kk[s][i] : q==2? bb[s][i] : res[s][i]);
        g->SetMarkerStyle(s?24:20); g->SetMarkerSize(0.9); g->SetMarkerColor(col[s]);
        g->SetLineColor(col[s]); g->SetLineWidth(2);
        if(!s){
          g->SetTitle(Form("%s   %s   %s",year,trig,ttl[q]));
          g->GetXaxis()->SetTitle("inst. lumi (/#mub/s)");
          g->GetYaxis()->SetTitle(""); g->GetXaxis()->SetLimits(0,0.0075);
          const double pad=0.12*(hi-lo)+1e-9;
          g->SetMinimum(lo-pad); g->SetMaximum(hi+pad);
          g->Draw("APL");
        } else g->Draw("PL same");
        lg->AddEntry(g,Form("ZDC%s",sn[s]),"pl");
      }
      lg->Draw();
      if(q==2){ TLine* z=new TLine(0,0,0.0075,0); z->SetLineStyle(2); z->SetLineColor(kGray+2); z->Draw(); }
    }
    c->SaveAs(Form("%s/lumi_dep/lumidep_%s_%s_%s.png",FIGD3,year,trig,nick[fig]));
    delete c;
  }
}
