// mu1(L): the 1n peak position against instantaneous luminosity, with fits.
// This is the only quantity an Xn (>=1 neutron) selection needs: a single
// multiplicative correction ZDC * mu_ref/mu1(L) fixes the peak, and the valley
// below it moves with the peak.  It is also the best-measured quantity here --
// one Gaussian on the strongest feature in the spectrum.
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
#include "TSystem.h"
#include "TString.h"
#include <vector>
#include <algorithm>

static const char* FIGD5="/Users/wl33/Library/CloudStorage/Dropbox/UPC/ZDC/figures";

void muFit(const char* year, const char* trig="ZDCOR",
           double Lmin=0.0010, double Lmax=0.0066, double nmin=20000){
  gStyle->SetOptStat(0); gStyle->SetTitleFontSize(0.050);
  gStyle->SetTitleAlign(23); gStyle->SetTitleX(0.56); gStyle->SetTitleY(0.995);
  gSystem->mkdir(Form("%s/mu1",FIGD5),kTRUE);
  TFile* f=TFile::Open(Form("lumispectra_%s.root",year));
  if(!f||f->IsZombie()){ printf("no spectra for %s\n",year); return; }
  TProfile* pL=(TProfile*)f->Get(Form("pL_%s_%s",year,trig));
  const char* sn[2]={"+","-"};
  TCanvas* c=new TCanvas(Form("cmu_%s",year),"",1500,660); c->Divide(2,1,0.004,0.004);
  printf("\n===== %s  %s   mu1(L), fit range %.4f - %.4f\n",year,trig,Lmin,Lmax);
  for(int s=0;s<2;s++){
    TH2D* h=(TH2D*)f->Get(Form("h%c_%s_%s",s?'M':'P',year,trig));
    if(!h) continue;
    std::vector<double> x,ex,mu,emu;
    double seed=0;
    { TH1D* a=h->ProjectionY("am");
      int bm=0; double v=0;
      for(int b=a->FindBin(1400.);b<=a->GetNbinsX();b++) if(a->GetBinContent(b)>v){v=a->GetBinContent(b);bm=b;}
      seed=bm? a->GetBinCenter(bm):2500; }
    for(int bx=1;bx<=h->GetNbinsX();bx++){
      TH1D* cc=h->ProjectionY(Form("m%s%d%d",year,s,bx),bx,bx);
      if(cc->Integral(cc->FindBin(1200.),cc->GetNbinsX())<nmin) continue;
      const double w=0.22*seed;
      int b1=cc->FindBin(seed-w), b2=cc->FindBin(seed+w), bm=0; double mx=0;
      for(int b=std::max(b1,cc->FindBin(1200.));b<=b2;b++) if(cc->GetBinContent(b)>mx){mx=cc->GetBinContent(b);bm=b;}
      if(!bm) continue;
      const double p=cc->GetBinCenter(bm);
      TF1 g("gm","gaus",p*0.74,p*1.30);
      g.SetParameters(mx,p,0.18*p); g.SetParLimits(1,p*0.82,p*1.18); g.SetParLimits(2,0.04*p,0.42*p);
      if(cc->Fit(&g,"QNR")!=0) continue;
      const double m=g.GetParameter(1), e=g.GetParError(1);
      if(m<500||m>40000||e<=0) continue;
      seed=m;
      const double Lv = pL->GetBinEntries(bx)>0? pL->GetBinContent(bx) : h->GetXaxis()->GetBinCenter(bx);
      if(Lv<Lmin||Lv>Lmax) continue;
      x.push_back(Lv); ex.push_back(0); mu.push_back(m); emu.push_back(e);
    }
    if(x.size()<5) continue;
    TPad* p=(TPad*)c->cd(s+1); p->SetLeftMargin(0.14); p->SetRightMargin(0.04); p->SetTopMargin(0.09);
    TGraphErrors* g=new TGraphErrors(x.size(),&x[0],&mu[0],&ex[0],&emu[0]);
    g->SetMarkerStyle(20); g->SetMarkerSize(0.8); g->SetMarkerColor(kBlack); g->SetLineColor(kBlack);
    g->SetTitle(Form("%s   1n peak #mu_{1}   ZDC%s",year,sn[s]));
    g->GetXaxis()->SetTitle("inst. lumi (/#mub/s)"); g->GetYaxis()->SetTitle("");
    g->GetXaxis()->SetLimits(0,0.0072); g->Draw("AP");
    TF1* f1=new TF1(Form("m1_%d",s),"pol1",Lmin,Lmax);
    TF1* f2=new TF1(Form("m2_%d",s),"pol2",Lmin,Lmax);
    TF1* f3=new TF1(Form("m3_%d",s),"pol3",Lmin,Lmax);
    g->Fit(f1,"QNR"); g->Fit(f2,"QNR"); g->Fit(f3,"QNR");
    f1->SetLineColor(kAzure+2); f1->SetLineStyle(2); f1->SetLineWidth(2); f1->Draw("same");
    f2->SetLineColor(kRed+1);   f2->SetLineWidth(2); f2->Draw("same");
    f3->SetLineColor(kGreen+3); f3->SetLineStyle(7); f3->SetLineWidth(2); f3->Draw("same");
    const double c1=f1->GetChisquare()/std::max(1,f1->GetNDF());
    const double c2=f2->GetChisquare()/std::max(1,f2->GetNDF());
    auto dev=[&](TF1* ff){ double d=0; for(size_t i=0;i<x.size();i++) d=std::max(d,fabs(mu[i]-ff->Eval(x[i]))/mu[i]); return 100*d; };
    TLegend* lg=new TLegend(0.18,0.70,0.44,0.88);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.036);
    lg->AddEntry(f1,"pol1","l");
    lg->AddEntry(f2,"pol2","l");
    lg->AddEntry(f3,"pol3","l");
    lg->Draw();
    // the largest deviation of each fit from the points, in percent -- what the
    // correction's residual would be if that form were used
    printf("%s ZDC%s  (%d pts, mu1 %.0f -> %.0f)\n",year,sn[s],(int)x.size(),mu.front(),mu.back());
    printf("   pol1  chi2/ndf %7.1f   max deviation %.2f%%   %9.1f %+10.4e*L\n",c1,dev(f1),f1->GetParameter(0),f1->GetParameter(1));
    printf("   pol2  chi2/ndf %7.1f   max deviation %.2f%%   %9.1f %+10.4e*L %+10.4e*L^2\n",c2,dev(f2),f2->GetParameter(0),f2->GetParameter(1),f2->GetParameter(2));
    printf("   pol3  chi2/ndf %7.1f   max deviation %.2f%%\n",f3->GetChisquare()/std::max(1,f3->GetNDF()),dev(f3));
    printf("   POL3 %s %s %s  %.6e %.6e %.6e %.6e  range %.5f %.5f  muRef(0.003) %.2f\n",
           year,trig,s?"M":"P",f3->GetParameter(0),f3->GetParameter(1),f3->GetParameter(2),f3->GetParameter(3),
           Lmin,Lmax,f3->Eval(0.003));
  }
  c->SaveAs(Form("%s/mu1/mu1_%s_%s.png",FIGD5,year,trig));
  delete c;
}
