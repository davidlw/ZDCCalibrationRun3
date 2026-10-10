// Skim two runs of one dataset into a single small file.  Two runs with
// different luminosity profiles are the control experiment: if the 1n peak and
// its width land on the same curve against L in both, then L is what sets them.

// Reuse a skim that already exists: these cost a full read of a multi-GB tree,
// and the same run is often wanted again for a different figure.
static bool skimExists(const char* out, const char* what){
  TFile* f=TFile::Open(out,"READ");
  if(f && !f->IsZombie()){
    TTree* t=(TTree*)f->Get("zdc_ana/zdcTree");
    if(t && t->GetEntries()>0){
      printf("%s already exists (%lld events) -- reusing, %s not re-skimmed\n",out,t->GetEntries(),what);
      f->Close(); return true;
    }
  }
  if(f) f->Close();
  return false;
}

void skim2(const char* inFile, const char* outFile, int runA, int runB){
  if(skimExists(outFile,Form("runs %d,%d",runA,runB))) return;
  ROOT::EnableImplicitMT(4);
  ROOT::RDataFrame d("zdc_ana/zdcTree", inFile);
  auto n = d.Filter(Form("run == %d || run == %d",runA,runB))
             .Snapshot("zdc_ana/zdcTree", outFile)->Count();
  printf("%s -> %s : %llu events\n",inFile,outFile,*n);
}
