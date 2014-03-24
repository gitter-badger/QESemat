#include <string.h>

//------------------------------------------------------
// some useful functions for fortran interface 
	char* fchar(const char* cchar, const int len){
		//function to create a string for fortran
		char *result=new char[len];
		memset(result,0,len); //fill string with spaces
		strcpy(result,cchar);
		return result;
	};


//------------------------------------------------------
// fortran functions to be used
extern "C"{
	int flux_init_();

	int flux_calc_spline_(int* NuAnu,int* Flavor);
	int flux_close_file_();
	int flux_has_table_(int* NuAnu,int* Flavor);
	int flux_open_file_(char *fname);
	int flux_print_table_(int* NuAnu,int* Flavor);
	int flux_read_hdr_();
	int flux_read_table_();

	double flux_get_df_(int* NuAnu,int* Flavor, double* Enu);
	double flux_get_emax_(int* NuAnu,int* Flavor);
	double flux_get_emin_(int* NuAnu,int* Flavor);
	double flux_get_zmax_(int* NuAnu,int* Flavor);
	double flux_get_zmin_(int* NuAnu,int* Flavor);
}
/// small interface
class Flux{

public:
	enum eNuAnu{kNu=1,kAnu=2};
	enum eFlavor{kE=1,kMu=2,kTau=3};
public:
	Flux():NuAnu(kNu),Flavor(kE){Init();}
	inline bool Init(){return flux_init_();}
	inline bool OpenFile(const char *fname){return flux_open_file_(fchar(fname,80));};
	inline bool ReadHead(){ return flux_read_hdr_();};
	inline bool ReadTable(){ return flux_read_table_();};
	inline bool CloseFile(){ return flux_close_file_();};
	inline bool CalcSpline(){ return flux_calc_spline_(&NuAnu, &Flavor);};

	inline bool HasTable(eNuAnu n,eFlavor f){ return flux_has_table_((int*)&n, (int*)&f);};
	inline bool PrintTable(eNuAnu n,eFlavor f){ return flux_print_table_((int*)&n, (int*)&f);};
	inline bool PrintTable(){ return flux_print_table_(&NuAnu, &Flavor);};

	double dF(double Enu){return flux_get_df_(&NuAnu,&Flavor, &Enu);};
	double Emax(){return flux_get_emax_( &NuAnu, &Flavor);};
	double Emin(){return flux_get_emin_( &NuAnu, &Flavor);};
	double Zmax(){return flux_get_zmax_( &NuAnu, &Flavor);};
	double Zmin(){return flux_get_zmin_( &NuAnu, &Flavor);};

	void SetNeutrino(eNuAnu n,eFlavor f){NuAnu=n; Flavor=f;}
public:
	int NuAnu, Flavor;
};

int main(int argc, char const *argv[])
{
	Flux f;
	f.OpenFile(argv[1]);
	f.ReadHead();
	f.ReadTable();
	f.SetNeutrino(Flux::kAnu,Flux::kE);
	f.PrintTable();
	return 0;
}
