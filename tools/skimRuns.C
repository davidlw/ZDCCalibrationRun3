// Skim the representative run of one dataset into a small file.
// Everything about the lumi-section figures -- binning, trigger split, which
// run -- can then be re-cut in seconds instead of re-reading a 12 GB tree.

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

void skimRuns(const char* inFile, const char* outFile, int runSel){
  if(skimExists(outFile,Form("run %d",runSel))) return;
  // Cap the threads: the skim is I/O bound, and several of these at once with
  // 8 workers each drove the load average to 60+ without going faster.
  ROOT::EnableImplicitMT(4);
  ROOT::RDataFrame d("zdc_ana/zdcTree", inFile);
  auto n = d.Filter(Form("run == %d",runSel)).Snapshot("zdc_ana/zdcTree", outFile)->Count();
  printf("%s -> %s : %llu events\n",inFile,outFile,*n);
}
