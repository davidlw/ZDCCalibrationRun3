// 1n peak per run, at fixed luminosity, to locate the ZDC calibration changes.
// Holding L fixed removes the luminosity drift, so what is left is calibration.
void runscale(const char* inFile, const char* tag, double L1=0.0015, double L2=0.0035){
  static std::map<ULong64_t,double> L; L.clear();
  { std::ifstream in("lumi_byls_2024HI.csv"); std::string line;
    while(std::getline(in,line)){ if(line.empty()||line[0]=='#') continue;
      std::vector<std::string> p; std::stringstream ss(line); std::string t;
      while(std::getline(ss,t,',')) p.push_back(t); if(p.size()<7) continue;
      size_t c0=p[0].find(':'), c1=p[1].find(':'); if(c0==std::string::npos||c1==std::string::npos) continue;
      double d=0; try{ d=std::stod(p[5]); }catch(...){ continue; }
      L[(ULong64_t)std::stoul(p[0].substr(0,c0))*100000ULL+std::stoul(p[1].substr(0,c1))]=d/23.31; } }
  ROOT::EnableImplicitMT(6);
  ROOT::RDataFrame d("zdc_ana/zdcTree",inFile);
  auto f = d.Filter("(trigBits & 0x8) != 0")
            .Define("Lv",[](UInt_t r,UInt_t l){ auto it=L.find((ULong64_t)r*100000ULL+l); return it==L.end()?-1.0:it->second; },{"run","lumi"})
            .Filter(Form("Lv > %f && Lv < %f",L1,L2));
  auto hm=f.Histo2D({"hm","",1100,387700,388800,240,0,16000},"run","zdcMinus");
  auto hp=f.Histo2D({"hp","",1100,387700,388800,240,0,16000},"run","zdcPlus");
  TH2D* HM=(TH2D*)hm.GetPtr(); TH2D* HP=(TH2D*)hp.GetPtr();
  printf("\n=== %s  (L %.4f-%.4f)\n%8s %9s %8s %8s\n",tag,L1,L2,"run","N","1n-","1n+");
  for(int b=1;b<=HM->GetNbinsX();b++){
    TH1D* c=HM->ProjectionY(Form("m%d",b),b,b);
    const double n=c->Integral(c->FindBin(1200.),c->GetNbinsX());
    if(n<20000) continue;
    auto pk=[](TH1* h){ int bm=0; double v=0;
      for(int i=h->FindBin(1400.);i<=h->GetNbinsX();i++) if(h->GetBinContent(i)>v){v=h->GetBinContent(i);bm=i;}
      return bm? h->GetBinCenter(bm):0.; };
    TH1D* cp=HP->ProjectionY(Form("p%d",b),b,b);
    printf("%8.0f %9.0f %8.0f %8.0f\n",HM->GetXaxis()->GetBinLowEdge(b),n,pk(c),pk(cp));
  }
}
