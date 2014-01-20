
void plot_file(char* fname, const int Ncomp=1,char* draw_opt=0){
  TString skp="";
  TString fmt="";
  TGraph* g[Ncomp];
  for(int N=0;N<Ncomp;++N){
    fmt="%lf ";
    fmt+=skp;
    fmt+="%lf";
    g[N]=new TGraph(fname,fmt);
    g[N]->SetTitle(Form("Component#%d",N));
    g[N]->SetLineWidth(2);
    g[N]->SetLineColor(N+1);
    if(draw_opt==0)g[N]->Draw((N>0)?"L":"AL");
    else g[N]->Draw(draw_opt);
    skp+="%*lf";
  }
    
}
