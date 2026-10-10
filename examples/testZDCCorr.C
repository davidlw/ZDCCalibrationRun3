// Closure test for ZDCLumiCorrection: the corrected 1n peak should be flat in
// luminosity.  Run from ~/ZDCtrees where the spectra files live.
#include "ZDCLumiCorrection.h"
void testZDCCorr(const char* year="2025", int runRep=400000, const char* csv="lumi_byls_2025HI.csv"){
  ZDCLumiCorrection c;
  c.addLumiTable(csv);
  TFile* f=TFile::Open(Form("lumispectra_%s.root",year));
  TProfile* pL=(TProfile*)f->Get(Form("pL_%s_ZDCOR",year));
  for(int s=0;s<2;s++){
    TH2D* h=(TH2D*)f->Get(Form("h%c_%s_ZDCOR",s?'M':'P',year));
    printf("\n%s ZDC%s   %8s %9s %9s %8s\n",year,s?"-":"+","L","mu1 raw","corrected","ratio");
    double first=0;
    for(int bx=1;bx<=h->GetNbinsX();bx+=3){
      TH1D* p=h->ProjectionY(Form("t%d%d",s,bx),bx,bx);
      if(p->Integral(p->FindBin(1200.),p->GetNbinsX())<20000) continue;
      int bm=0; double v=0;
      for(int b=p->FindBin(1400.);b<=p->GetNbinsX();b++) if(p->GetBinContent(b)>v){v=p->GetBinContent(b);bm=b;}
      const double mu=p->GetBinCenter(bm);
      const double L=pL->GetBinContent(bx);
      const double m1=c.mu1(runRep,(ZDCLumiCorrection::Side)s,L);
      const double cor=mu*c.muRef(runRep,(ZDCLumiCorrection::Side)s)/m1;
      if(!first) first=cor;
      printf("        %8.4f %9.0f %9.0f %8.3f\n",L,mu,cor,cor/first);
    }
  }
}
