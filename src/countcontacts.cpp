#include <mechsys/dem/domain.h>
#include <stdlib.h>
using namespace std;
using namespace DEM;
int main(int argc, char **argv) try
{
	string datafile="countcontacts.drp";
	string domainin;
	string domainout;
	ifstream datain;
	double Kn;
	double Kt;
	double Gn;
	double Gt;
	double Mu;
	double dt;
	double tf;
	double dtout;
	double Alpha;
	double Kinematicenergy;
	string filekey;
	int visualization;
	int processornumber;
	datain.open(datafile.c_str());
	datain >> domainin; 				datain.ignore(200,'\n');
	datain >> domainout;				datain.ignore(200,'\n');
	datain.close();
	DEM::Domain particles;
	particles.Load(domainin.c_str());
	datain.open("countcontacts.par");
	datain >> Kn;				datain.ignore(200,'\n');
	datain >> Kt;				datain.ignore(200,'\n');
	datain >> Gn;				datain.ignore(200,'\n');
	datain >> Gt;				datain.ignore(200,'\n');
	datain >> Mu;				datain.ignore(200,'\n');
	datain >> dt;				datain.ignore(200,'\n');
	datain >> tf;				datain.ignore(200,'\n');
	datain >> dtout;			datain.ignore(200,'\n');
	datain >> filekey;			datain.ignore(200,'\n');
	datain >> visualization;	datain.ignore(200,'\n');
	datain >> processornumber;	datain.ignore(200,'\n');
    particles.Solve  (/*tf*/tf, /*dt*/dt, /*dtOut*/dtout, NULL, NULL, /*filekey*/filekey.c_str(),/*Visit visualization*/visualization,/*N_proc*/processornumber, /*kinematic energy*/Kinematicenergy);
	int numbercontacts[1000];
	int surface[1000];
	for (int i=0; i<1000;i++)
	{
		numbercontacts[i]=0;
	}
	for (size_t i=0; i<particles.Particles.Size();i++)
	{
		numbercontacts[size_t(particles.Particles[i]->Cn)]++;
	}
	double average=0.0;
	int particle =0;
	for (int i=0; i<100;i++)
	{
		average += i*numbercontacts[i];
		particle +=numbercontacts[i];
	}
	average /= particle;
	ofstream dataout;																							// export to text file
	dataout.open(domainout.c_str());
	dataout << average << " Average number of contacts \n";
	dataout << particle << " Number of particles \n";
	for (int i =0; i<100; i++)
	{
		dataout << i<<" "<<numbercontacts[i] << "\n";
	}
	dataout << "\n";
	dataout.close();
}
MECHSYS_CATCH
