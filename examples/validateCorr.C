// Validation: apply ZDCLumiCorrection event by event and overlay the spectra
// from lumi-section slices of one run.  Raw, the five slices march apart as the
// luminosity decays; corrected, they should sit on top of each other.
#include "ZDCLumiCorrection.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>

void validateCorr(const char* skim="runskim/skim_2026_run404886_404730.root",
                  int run=404886, const char* csv="lumi_byls_2026HI.csv",
                  const char* year="2026", ULong64_t mask=0x8, int NSL=5)
{
  ZDCLumiCorrection corr;
  corr.addLumiTable(csv);
  TFile* f=TFile::Open(skim); TTree* t=(TTree*)f->Get("zdc_ana/zdcTree");

  // equal-statistics lumi-section slices
  std::vector<int> v;
  { TTreeReader r0(t);
    TTreeReaderValue<UInt_t> rn(r0,"run"), ls(r0,"lumi"); TTreeReaderValue<ULong64_t> tb(r0,"trigBits");
    while(r0.Next()) if((int)*rn==run && (*tb&mask)) v.push_back(*ls); }
  std::sort(v.begin(),v.end());
  std::vector<int> lo(NSL),hi(NSL); std::vector<double> mL(NSL,0); std::vector<long> nL(NSL,0);
  for(int k=0;k<NSL;k++){
    lo[k]=v[(size_t)(1.0*k/NSL*(v.size()-1))];
    hi[k]=(k==NSL-1)? v.back() : v[(size_t)(1.0*(k+1)/NSL*(v.size()-1))];
  }
  TH1D* h[2][2][8];                       // [side][raw/corr][slice]
  for(int s=0;s<2;s++) for(int c=0;c<2;c++) for(int k=0;k<NSL;k++)
    h[s][c][k]=new TH1D(Form("v%d%d%d",s,c,k),"",200,0,16000);

  TTreeReader r(t);
  TTreeReaderValue<UInt_t> rn(r,"run"), ls(r,"lumi");
  TTreeReaderValue<ULong64_t> tb(r,"trigBits");
  TTreeReaderValue<Float_t> zp(r,"zdcPlus"), zm(r,"zdcMinus");
  long nbad=0, ntot=0;
  while(r.Next()){
    if((int)*rn!=run || !(*tb&mask)) continue;
    int k=-1; for(int j=0;j<NSL;j++) if((int)*ls>=lo[j] && (int)*ls<=hi[j]){ k=j; break; }
    if(k<0) continue;
    ntot++;
    const double L=corr.instLumi(*rn,*ls);
    if(L>0){ mL[k]+=L; nL[k]++; }
    h[0][0][k]->Fill(*zp); h[1][0][k]->Fill(*zm);
    const double cp=corr.correct(*rn,*ls,ZDCLumiCorrection::Plus,*zp);
    if(!corr.valid()) nbad++;
    const double cm=corr.correct(*rn,*ls,ZDCLumiCorrection::Minus,*zm);
    h[0][1][k]->Fill(cp); h[1][1][k]->Fill(cm);
  }
  for(int k=0;k<NSL;k++) if(nL[k]) mL[k]/=nL[k];
  printf("%ld events, %ld (%.2f%%) with no luminosity -- left uncorrected\n",ntot,nbad,100.*nbad/ntot);

  gStyle->SetOptStat(0); gStyle->SetTitleFontSize(0.052);
  gStyle->SetTitleAlign(23); gStyle->SetTitleX(0.56); gStyle->SetTitleY(0.995);
  gSystem->mkdir("/Users/wl33/Library/CloudStorage/Dropbox/UPC/ZDC/figures/validation",kTRUE);
  const int col[6]={kAzure+2,kTeal+3,kSpring-6,kOrange+7,kRed+1,kMagenta+2};
  TCanvas* c=new TCanvas("cv","",1500,900); c->Divide(2,2,0.004,0.004);
  for(int s=0;s<2;s++) for(int q=0;q<2;q++){
    TPad* p=(TPad*)c->cd(1+q+2*s); p->SetLeftMargin(0.11); p->SetRightMargin(0.04); p->SetTopMargin(0.085);
    double mx=0, xhi=0;
    for(int k=0;k<NSL;k++){ TH1D* hh=h[s][q][k];
      if(hh->Integral()>0) hh->Scale(1./hh->Integral()); }
    // peak of the first slice sets the axis
    { TH1D* hh=h[s][q][0]; int bm=0; double vv=0;
      for(int b=hh->FindBin(1400.);b<=hh->GetNbinsX();b++) if(hh->GetBinContent(b)>vv){vv=hh->GetBinContent(b);bm=b;}
      xhi=3.6*hh->GetBinCenter(bm); }
    for(int k=0;k<NSL;k++) for(int b=h[s][q][k]->FindBin(1500.);b<=h[s][q][k]->FindBin(xhi);b++)
      mx=std::max(mx,h[s][q][k]->GetBinContent(b));
    TLegend* lg=new TLegend(0.46,0.62,0.97,0.88);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.031);
    for(int k=0;k<NSL;k++){
      TH1D* hh=h[s][q][k];
      hh->SetMarkerStyle(20); hh->SetMarkerSize(0.5); hh->SetMarkerColor(col[k%6]);
      hh->SetLineColor(col[k%6]); hh->SetLineWidth(1);
      hh->GetXaxis()->SetRangeUser(600.,xhi); hh->GetYaxis()->SetRangeUser(0,mx*1.6);
      hh->SetTitle(Form("%s (run %d)   ZDC%s   %s",year,run,s?"-":"+",q?"corrected":"raw"));
      hh->GetXaxis()->SetTitle(Form("ZDC%s (a.u.)",s?"-":"+"));
      hh->GetYaxis()->SetTitle("");
      hh->Draw(k?"P same":"P");
      int bm=0; double vv=0;
      for(int b=hh->FindBin(1400.);b<=hh->GetNbinsX();b++) if(hh->GetBinContent(b)>vv){vv=hh->GetBinContent(b);bm=b;}
      const double pk=hh->GetBinCenter(bm);
      TF1* g=new TF1(Form("gv%d%d%d",s,q,k),"gaus",pk*0.74,pk*1.30);
      g->SetParameters(vv,pk,0.2*pk); hh->Fit(g,"QNR");
      g->SetLineColor(kBlack); g->SetLineWidth(2); g->Draw("same");
      lg->AddEntry(hh,Form("L = %.4f,  1n = %.0f",mL[k],g->GetParameter(1)),"p");
    }
    lg->Draw();
  }
  c->SaveAs(Form("/Users/wl33/Library/CloudStorage/Dropbox/UPC/ZDC/figures/validation/corr_%s_run%d.png",year,run));
  delete c;
}
