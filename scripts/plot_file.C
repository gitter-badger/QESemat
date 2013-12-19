
void plot_file(char* fname){
  TGraph *g0=new TGraph(fname,"%lf %lf");
  g0->SetTitle("Total"); g0->SetLineWidth(2); g0->SetLineColor(kBlack);
  TGraph *g1=new TGraph(fname,"%lf %*lf %lf");
  g1->SetTitle("Component#1"); g1->SetLineWidth(2); g1->SetLineColor(2);
  TGraph *g2=new TGraph(fname,"%lf %*lf %*lf %lf");
  g2->SetTitle("Component#2"); g2->SetLineWidth(2); g2->SetLineColor(3);
  TGraph *g3=new TGraph(fname,"%lf %*lf %*lf %*lf %lf");
  g3->SetTitle("Component#3"); g3->SetLineWidth(2); g3->SetLineColor(4);
  TGraph *g4=new TGraph(fname,"%lf %*lf %*lf %*lf %*lf %lf");
  g4->SetTitle("Component#4"); g4->SetLineWidth(2); g4->SetLineColor(5);
  TGraph *g5=new TGraph(fname,"%lf %*lf %*lf %*lf %*lf %*lf %lf");
  g5->SetTitle("Component#5"); g5->SetLineWidth(2); g5->SetLineColor(6);
  g0->Draw("AL");
  g1->Draw("L");
  g2->Draw("L");
  g3->Draw("L");
  g4->Draw("L");
  g5->Draw("L");


  
}
