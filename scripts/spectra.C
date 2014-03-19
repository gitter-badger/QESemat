#include "TStyle.h"

void SetStyle(){
	gStyle->SetFillColor(0);
	gStyle->SetFrameLineWidth(1);
	gStyle->SetPadLeftMargin(0.075);
	gStyle->SetPadRightMargin(0.007);
	gStyle->SetPadTopMargin(0.005);		
	gStyle->SetPadBottomMargin(0.08);
	gStyle->SetTitleOffset(1.05,"Y");
	gStyle->SetTitleOffset(1.15,"X");
//	gStyle->SetTextFont(132);
}

void spectra(){
	
	const int chf=2,cht=2,che=770,cha=2;
	const char *nid="../../input/",*nod="../../resnlog/";
	const char nsp[10]="nova_flux",f[chf]={'e','m'},t[cht]={'n','a'};
	const char fn[chf][4]={"e","#mu"},tn[cht][8]={"#nu","#bar#nu"};
	char nof[256],pfn[256];
	FILE *f140=0;
	double sp[chf][cht][che][cha],x,y;
	const int col[chf][cht]={{4,6},{1,2}},sh[chf][cht]={{9,3},{1,7}},w[chf][cht]={{1,2},{1,2}};
	int e,a,nf,nt,n,ndr;
	
	TCanvas c("c","Canvas",1200,800);
	TGraph grsp[chf][cht];
	
	SetStyle();

    for(nf=0; nf<chf; nf++){
        for(nt=0; nt<cht; nt++){
            sprintf(pfn,"%s_%c%c",nsp,f[nf],t[nt]);
            sprintf(nof,"%s%s.sng",nid,pfn);
            f140=fopen(nof,"r");
            for(e=0; e<che; e++){
                for(a=0; a<cha; a++){
                    fscanf(f140,"%le",&(sp[nf][nt][e][a]));
                }
            }
            fclose(f140);
        }
    }
	
	for(ndr=0; ndr<2; ndr++){
        c.cd();
        for(nf=0; nf<chf; nf++){
            for(nt=0; nt<cht; nt++){
                for(e=0; e<che; e++){
                    x=sp[nf][nt][e][0];
                    y=sp[nf][nt][e][1];
                    grsp[nf][nt].SetPoint(e,x,y);
                }
            }
        }
        grsp[1][0].GetXaxis()->SetTitle("#font[132]{E_{#nu} (GeV)}");
        grsp[1][0].GetXaxis()->SetRangeUser(0,40);
        grsp[1][0].GetYaxis()->SetTitle("#font[132]{dF_{#nu}/dE_{#nu} (???)}");//GeV^{-1} s^{-1} cm^{-2}
        grsp[1][0].GetYaxis()->SetRangeUser(1e-1,2e5);
        grsp[1][0].Draw("AC");
        for(nf=0; nf<chf; nf++){
            for(nt=0; nt<cht; nt++){
                grsp[nf][nt].SetLineColor(col[nf][nt]);
                grsp[nf][nt].SetLineStyle(sh[nf][nt]);
                grsp[nf][nt].SetLineWidth(w[nf][nt]);
                grsp[nf][nt].Draw("CSAME");
            }
        }
        nt=0;
        TLegend* legl=new TLegend(0.24,0.82,0.47,0.95);
        for(nf=0; nf<chf; nf++){
            legl->AddEntry(&grsp[nf][nt],Form("#font[132]{%s_{%s}}",tn[nt],fn[nf]),"L");
        }
        nt=1;
        TLegend* legr=new TLegend(0.33,0.82,0.56,0.95);
        for(nf=0; nf<chf; nf++){
            legr->AddEntry(&grsp[nf][nt],Form("#font[132]{%s_{%s}}",tn[nt],fn[nf]),"L");
        }
        legl->SetBorderSize(0);
        legl->Draw();
        legr->SetBorderSize(0);
        legr->Draw();
        c.SetLogy();
        sprintf(pfn,"%s",nsp);
        sprintf(nof,"%s%s.eps",nod,pfn);
        c.Print((const char*)nof);
        delete legl;
        delete legr;
    }
}
