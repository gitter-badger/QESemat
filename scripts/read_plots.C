void SetStyle(){
	gStyle->SetLineWidth(2);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetTitleYOffset(1.3);
}

TGraph* ReadFile(char* fname, int ncomp, char* gname, char* gtitle=0)
{
	if(gSystem->AccessPathName(fname))return 0;
	if(gtitle==0)gtitle=gname;
	TString fmt="%lf";
	for (int i = 0; i < ncomp; ++i){
		fmt+=(i==ncomp-1)?" %lf":" %*lf";
	}
	// printf("fmt[%d]=%s\n",ncomp,fmt.Data());
	TGraph* graph=new TGraph(fname,fmt);
	graph->SetName(gname);
	return graph;
	//divide everything by Enu
	double x,y;
	for(unsigned j = 0; j < graph->GetN(); ++j) {
		/* code */
		graph->GetPoint(j,x,y);
		graph->SetPoint(j,x,y/x);
	}
	return graph;
}

TString nutype;
TString Description="";

TString gen_title(char* fname){
	TString s=fname;
	int idx;
	idx=s.Last('_');
	nutype=s(idx+1,2);
	printf("nutype=%s\n",nutype.Data());
	TString title;
	title=(nutype(0)=='n')?"#nu":"#bar{#nu}";

	switch (nutype(1)){
		case 'e': title+="_{e}"; break;
		case 'm': title+="_{#mu}"; break;
		case 't': title+="_{#tau}"; break;
	}
	title+=Description;
	// char corv=s(idx+1);
	// title+=", M_{A}";
	// switch (corv){
	// 	case 'c':title+=" constant"; break;
	// 	case 'v':title+=" effective"; break;
	// }
	if(s.Contains("flux"))title+="; P_{l}, GeV; N_{events} per kg ";
	else title+="; E_{#nu}, GeV; #sigma_{QES}#times 10^{38}";
	return title;
}

int Ncomp=1; //component to read

void read_files(char* file_base="Scintillator_ae", int color=kBlue){
	TGraph* g=0;
	// setup all suffixes here
	char* suffix[]={"0.95","1.1","1.2","1.35","0.95"};
	char corv[]  ="ccccv";
	int  styles[]={1,2,3,7,1};
	int colors[]={color,color,color,color,kBlack};
	// done
	unsigned Nsuf=sizeof(suffix)/sizeof(char*);
	TMultiGraph* mg=new TMultiGraph(file_base,gen_title(file_base));
	TLegend* leg=new TLegend(0.3,0.2,0.5,0.4);

	int Ngr=0; //number of graphs

	for(unsigned i = 0; i < Nsuf; ++i) {
		char* fname=Form("%s_%c_%s.dat",file_base,corv[i],suffix[i]);
		printf("read file %s\n",fname);
		g=ReadFile(fname,Ncomp,"test");
		if(!g)continue;
		g->SetLineStyle(styles[i]);
		g->SetLineColor(colors[i]);
		mg->Add(g);
		if(corv[i]=='c')leg->AddEntry(g,Form("M_{A}=%s",suffix[i]),"l");
		else{
			g->SetLineWidth(3);
			leg->AddEntry(g,"M_{A} effective","l");
		}
		Ngr++;
	}
	mg->SetMinimum(0.0);
	mg->Draw("al");
	if(Ngr>1)leg->Draw();
}

void read_plots(char* fbase="Scintillator", int ncomp=1, char* desc=""){
	SetStyle();
	Description=desc;
	Ncomp=ncomp;
	TString s(fbase);
	if(s.Contains("flux"))gStyle->SetOptLogy(1);
	else gStyle->SetOptLogx(1);
	TCanvas* c=0;
	c=new TCanvas("E","E",800,400);
	c->Divide(2,1);
	c->cd(1); read_files(Form("%s_ne",fbase),kRed+2);
	c->cd(2); read_files(Form("%s_ae",fbase),kRed+2);
	c->Print(Form("%s_E_n%d.eps",fbase,ncomp));

 	c=new TCanvas("mu","mu",800,400);
	c->Divide(2,1);
	c->cd(1); read_files(Form("%s_nm",fbase),kBlue+2);
	c->cd(2); read_files(Form("%s_am",fbase),kBlue+2);
	c->Print(Form("%s_M_n%d.eps",fbase,ncomp));	
}
