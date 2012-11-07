#include <mechsys/dem/domain.h>
using namespace std;
using namespace DEM;
int main(int argc, char **argv) try
{
	int count1 =0;
	int count2 =0;
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
	double roundratio;
	DEM::Domain particles;
	//DEM::Domain particlesbefore;
	Array<Vec3_t> posbefore;
	string datafile="compact.drp";
	string domainin;
	string domainout;
	int tag;
	String axis;
	double force;
	ifstream datain;
	datain.open(datafile.c_str());
	datain >> domainin; 				datain.ignore(200,'\n');
	datain >> domainout;				datain.ignore(200,'\n');
	datain >> tag;						datain.ignore(200,'\n');
	datain >> axis;						datain.ignore(200,'\n');
	datain >> force;					datain.ignore(200,'\n');
	datain.close();
	particles.Load(domainin.c_str());
	datain.open("compact.par");
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
	datain >> roundratio;		datain.ignore(200,'\n');
	datain >> Alpha;     		datain.ignore(200,'\n');
	datain >> Kinematicenergy;  datain.ignore(200,'\n');
	datain.close();
	cout << "data reading \n";
    for (size_t i=-1;i>-7;i--)
    {
        particles.GetParticle(i)->FixVeloc();
    }
	if (axis=="z")
	{
		particles.GetParticle(tag)->vzf = false;
		particles.GetParticle(tag)->Ff=Vec3_t(0,0,force);
	}
    Dict P;
    for (size_t i=0;i<20;i++)
    {
        P.Set(/*Tag*/i,"Kn Kt Gn Gt Mu",Kn,Kt,Gn,Gt,Mu);
    }
    cout << "set parameters \n";
    Vec3_t g(0.0,0.0,-9.8);
    particles.Alpha = Alpha;
    for (size_t i=0;i<particles.Particles.Size();i++)
        {
            particles.Particles[i]->Ff = particles.Particles[i]->Props.m*g;
            posbefore.Push(particles.Particles[i]->x);
        }
    particles.Solve  (/*tf*/tf, /*dt*/dt, /*dtOut*/dtout, NULL, NULL, /*filekey*/filekey.c_str(),/*Visit visualization*/visualization,/*N_proc*/processornumber, /*kinematic energy*/Kinematicenergy);
    particles.Save (domainout.c_str());
    cout << "solve domain \n";
	return 0;
}
MECHSYS_CATCH
