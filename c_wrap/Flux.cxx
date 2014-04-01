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
	int Flux_Init();

	int Flux_Calc_Spline(int* NuAnu,int* Flavor);
	int Flux_Close_File();
	int Flux_Has_Table(int* NuAnu,int* Flavor);
	int Flux_Open_File(char *fname);
	int Flux_Print_Table(int* NuAnu,int* Flavor);
	int Flux_Read_Hdr();
	int Flux_Read_Table();

	double Flux_Get_dF(int* NuAnu,int* Flavor, double* Enu);
	double Flux_Get_Emax(int* NuAnu,int* Flavor);
	double Flux_Get_Emin(int* NuAnu,int* Flavor);
	double Flux_Get_Zmax(int* NuAnu,int* Flavor);
	double Flux_Get_Zmin(int* NuAnu,int* Flavor);
}
/// small interface
class Flux{

public:
	enum eNuAnu{kNu=1,kAnu=2};
	enum eFlavor{kE=1,kMu=2,kTau=3};
public:
	Flux():NuAnu(kNu),Flavor(kE){Init();}
	inline bool Init(){return Flux_Init();}
	inline bool OpenFile(const char *fname){return Flux_Open_File(fchar(fname,80));};
	inline bool ReadHead(){ return Flux_Read_Hdr();};
	inline bool ReadTable(){ return Flux_Read_Table();};
	inline bool CloseFile(){ return Flux_Close_File();};
	inline bool CalcSpline(){ return Flux_Calc_Spline(&NuAnu, &Flavor);};

	inline bool HasTable(eNuAnu n,eFlavor f){ return Flux_Has_Table((int*)&n, (int*)&f);};
	inline bool PrintTable(eNuAnu n,eFlavor f){ return Flux_Print_Table((int*)&n, (int*)&f);};
	inline bool PrintTable(){ return Flux_Print_Table(&NuAnu, &Flavor);};

	double dF(double Enu){return Flux_Get_dF(&NuAnu,&Flavor, &Enu);};
	double Emax(){return Flux_Get_Emax( &NuAnu, &Flavor);};
	double Emin(){return Flux_Get_Emin( &NuAnu, &Flavor);};
	double Zmax(){return Flux_Get_Zmax( &NuAnu, &Flavor);};
	double Zmin(){return Flux_Get_Zmin( &NuAnu, &Flavor);};

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
