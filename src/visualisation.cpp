#include <mechsys/dem/domain.h>
#include <stdlib.h>
using namespace std;
using namespace DEM;
int main(int argc, char **argv) try
{
	string datafile="visualisation.drp";
	string domainin;
	string domainout;
	ifstream datain;
	datain.open(datafile.c_str());
	datain >> domainin; 				datain.ignore(200,'\n');
	datain >> domainout;				datain.ignore(200,'\n');
	datain.close();
	DEM::Domain particles;
	particles.Load(domainin.c_str());
	particles.Save(domainout.c_str());
	particles.WriteXDMF(domainout.c_str());// export to draw visual results
	cout << "Write: " <<domainin <<" to: " << domainout<<"\n";
}
MECHSYS_CATCH
