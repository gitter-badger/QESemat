
char* get_names(char* fname, int n){
  FILE* f=fopen(fname,"r");
  char line[256];
  char* tokens=new char[3*n];
  fgets(line,256,f);
  char* token=0;
  int idx=13;
  for (int i = 0; i < n; ++i)
  {
    /* code */
    idx+=21;
    if(idx>strlen(line))break;
    token=tokens+3*i;
    printf("line[%d-%d]=%c%c\n",idx,idx+1,line[idx],line[idx+1]);
    token[0]=line[idx];
    token[1]=line[idx+1];
    token[2]='\0';
    //strncpy(token,line+idx,2);
    printf("token=%s\n",token);
    idx+=2;
  }
  fclose(f);
  return tokens;
};

void plot_file(char* fname, const int Ncomp=1,char* suf="_",char* draw_opt="al"){
  if(Ncomp>1)char* tokens=get_names(fname,Ncomp-1);
  TString skp="";
  TString fmt="";
  TGraph* g[Ncomp];
  TMultiGraph* mgr=new TMultiGraph;
  mgr->SetTitle("Event rate; P_{lep}, GeV; N_{events}");
  TLegend * leg=new TLegend(0.8,0.8,1.,1.);
  TFile fout("Hist.root","UPDATE");
  fout.cd();
  for(int N=0;N<Ncomp;++N){
    fmt="%lf ";
    fmt+=skp;
    fmt+="%lf";
    g[N]=new TGraph(fname,fmt);
    if(N!=Ncomp-1)
      g[N]->SetTitle(tokens+3*N);
    else 
      g[N]->SetTitle("Total");
    g[N]->SetLineWidth(2);
    g[N]->SetLineColor(N+1);
    leg->AddEntry(g[N],0,"l");
    //if(draw_opt==0)g[N]->Draw((N>0)?"L":"AL");
    //else g[N]->Draw(draw_opt);
    g[N]->Write(Form("comp%s_%d",suf,N));
    mgr->Add(g[N]);
    skp+="%*lf";
  }
  fout.Close();
  mgr->Draw(draw_opt);
  leg->Draw();
}
