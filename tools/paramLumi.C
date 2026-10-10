// Parametrising b(L) and k(L), one year at a time.
//
// Points carry the statistical error from the Gaussian peak fits, propagated:
//   k = mu2 - mu1      -> err = sqrt(e1^2 + e2^2)
//   b = 2*mu1 - mu2    -> err = sqrt(4*e1^2 + e2^2)
// Without those the chi2 means nothing and every form looks equally good.
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include "TString.h"
#include <vector>
#include <algorithm>

static const char* FIGD4="/Users/wl33/Library/CloudStorage/Dropbox/UPC/ZDC/figures";

static bool pk1(TH1* h, double seed, double win, double lo, double hi,
                double& mu, double& e, double& sg){
  const double w=win*seed;
  int b1=h->FindBin(seed-w), b2=h->FindBin(seed+w), bm=0; double mx=0;
  for(int b=std::max(b1,1);b<=std::min(b2,h->GetNbinsX());b++)
    if(h->GetBinContent(b)>mx){ mx=h->GetBinContent(b); bm=b; }
  if(!bm) return false;
  const double p=h->GetBinCenter(bm);
  TF1 g("gq","gaus",p*lo,p*hi);
  g.SetParameters(mx,p,0.18*p);
  g.SetParLimits(1,p*(1-0.8*win),p*(1+0.8*win));
  g.SetParLimits(2,0.04*p,0.42*p);
  if(h->Fit(&g,"QNR")!=0) return false;
  mu=g.GetParameter(1); e=g.GetParError(1); sg=fabs(g.GetParameter(2));
  return mu>500 && mu<40000 && sg/mu<0.5 && e>0;
}

void paramLumi(const char* year, const char* trig="ZDCOR",
               double Lmin=0.0010, double Lmax=0.0066, double nmin=20000){
  gStyle->SetOptStat(0); gStyle->SetTitleFontSize(0.052);
  gStyle->SetTitleAlign(23); gStyle->SetTitleX(0.56); gStyle->SetTitleY(0.995);
  gSystem->mkdir(Form("%s/param",FIGD4),kTRUE);
  TFile* f=TFile::Open(Form("lumispectra_%s.root",year));
  if(!f||f->IsZombie()){ printf("no spectra for %s\n",year); return; }
  TProfile* pL=(TProfile*)f->Get(Form("pL_%s_%s",year,trig));
  const char* sn[2]={"+","-"};
  std::vector<double> x[2],ex[2],kk[2],ek[2],bb[2],eb[2];
  for(int s=0;s<2;s++){
    TH2D* h=(TH2D*)f->Get(Form("h%c_%s_%s",s?'M':'P',year,trig));
    if(!h){ printf("missing %s %s\n",year,trig); return; }
    double seed=0;
    { TH1D* a=h->ProjectionY("aq");
      int bm=0; double v=0;
      for(int b=a->FindBin(1400.);b<=a->GetNbinsX();b++) if(a->GetBinContent(b)>v){v=a->GetBinContent(b);bm=b;}
      seed=bm? a->GetBinCenter(bm):2500; }
    for(int bx=1;bx<=h->GetNbinsX();bx++){
      TH1D* c=h->ProjectionY(Form("q%s%d%d",year,s,bx),bx,bx);
      if(c->Integral(c->FindBin(1200.),c->GetNbinsX())<nmin) continue;
      double m1,e1,s1,m2,e2,s2;
      if(!pk1(c,seed,0.22,0.74,1.30,m1,e1,s1)) continue;
      if(!pk1(c,2*m1,0.16,0.87,1.13,m2,e2,s2)) continue;
      if(m2<1.5*m1||m2>2.5*m1) continue;
      seed=m1;
      const double Lv = pL->GetBinEntries(bx)>0? pL->GetBinContent(bx) : h->GetXaxis()->GetBinCenter(bx);
      if(Lv<Lmin||Lv>Lmax) continue;
      x[s].push_back(Lv); ex[s].push_back(0);
      kk[s].push_back(m2-m1);      ek[s].push_back(sqrt(e1*e1+e2*e2));
      bb[s].push_back(2*m1-m2);    eb[s].push_back(sqrt(4*e1*e1+e2*e2));
    }
  }
  TCanvas* c=new TCanvas(Form("cpar_%s",year),"",1500,900); c->Divide(2,2,0.004,0.004);
  printf("\n===== %s  %s   fit range L = %.4f - %.4f\n",year,trig,Lmin,Lmax);
  for(int q=0;q<2;q++) for(int s=0;s<2;s++){
    auto& y  = q? bb[s] : kk[s];
    auto& ey = q? eb[s] : ek[s];
    if(x[s].size()<5) continue;
    TPad* p=(TPad*)c->cd(1+s+2*q); p->SetLeftMargin(0.14); p->SetRightMargin(0.04); p->SetTopMargin(0.085);
    TGraphErrors* g=new TGraphErrors(x[s].size(),&x[s][0],&y[0],&ex[s][0],&ey[0]);
    g->SetMarkerStyle(20); g->SetMarkerSize(0.8); g->SetMarkerColor(kBlack); g->SetLineColor(kBlack);
    g->SetTitle(Form("%s   %s   ZDC%s",year,q?"offset b":"per-neutron scale k",sn[s]));
    g->GetXaxis()->SetTitle("inst. lumi (/#mub/s)"); g->GetYaxis()->SetTitle("");
    g->GetXaxis()->SetLimits(0,0.0072);
    g->Draw("AP");
    TF1* f1=new TF1(Form("f1_%d%d",q,s),"pol1",Lmin,Lmax);
    TF1* f2=new TF1(Form("f2_%d%d",q,s),"pol2",Lmin,Lmax);
    g->Fit(f1,"QNR"); g->Fit(f2,"QNR");
    f1->SetLineColor(kAzure+2); f1->SetLineWidth(2); f1->SetLineStyle(2); f1->Draw("same");
    f2->SetLineColor(kRed+1);   f2->SetLineWidth(2); f2->Draw("same");
    const double c1=f1->GetChisquare()/std::max(1,f1->GetNDF());
    const double c2=f2->GetChisquare()/std::max(1,f2->GetNDF());
    TLegend* lg=new TLegend(0.16,0.72,0.72,0.88);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.036);
    lg->AddEntry(f1,Form("linear   #chi^{2}/ndf = %.1f",c1),"l");
    lg->AddEntry(f2,Form("quadratic   #chi^{2}/ndf = %.1f",c2),"l");
    lg->Draw();
    if(q) { TLine* z=new TLine(0,0,0.0072,0); z->SetLineStyle(3); z->SetLineColor(kGray+2); z->Draw(); }
    printf("%s ZDC%s  %s\n",year,sn[s],q?"b":"k");
    printf("   linear    : %9.1f %+10.4e*L                  chi2/ndf %8.1f  (%d pts)\n",
           f1->GetParameter(0),f1->GetParameter(1),c1,(int)x[s].size());
    printf("   quadratic : %9.1f %+10.4e*L %+10.4e*L^2   chi2/ndf %8.1f\n",
           f2->GetParameter(0),f2->GetParameter(1),f2->GetParameter(2),c2);
  }
  c->SaveAs(Form("%s/param/param_%s_%s.png",FIGD4,year,trig));
  delete c;
}
