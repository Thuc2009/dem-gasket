#include <mechsys/dem/domain.h>
using namespace std;
using namespace DEM;
void Report (DEM::Domain & particles, void *UD)
{
	String fn;
	fn.Printf    ("%s_%s_%04d", particles.FileKey.CStr(), "bf",particles.idx_out);
	particles.Save (fn.c_str());
}
int main(int argc, char **argv) try
{
	string filename;
	string filename1;
	if (argc<2)
		{
			filename = "continuousdrop.drp";
			filename1 = "continuousdrop.par";
		}
	else if (argc <3)
		{
			filename = argv[1];
			filename1 = "continuousdrop.par";
		}
	else
		{
			filename =argv[1];
			filename1 =argv[2];
		}
	int count1 =0;
	int count2 =0;
	double Kn;
	double Kt;
	double Gn;
	double Gt;
	double Gv;
	double Mu;
	double dt;
	double tf;
	double dtout;
	double Alpha;
	double Kinematicenergy;
	double roundratio;

	int starttag;
	int numbershapes;
	int numberintervals;
	string filekey;
	int visualization;
	int processornumber;

	DEM::Domain particles;
	Vec3_t Xmin;
	Vec3_t Xmax;
	double size[3];
	string domainin;
	string domainout;
	string boundary;
	ifstream datain;
	datain.open(filename.c_str());
	datain >> domainin; 				datain.ignore(200,'\n');
	datain >> domainout;				datain.ignore(200,'\n');
	datain >> boundary;					datain.ignore(200,'\n');
	datain.close();
	particles.Load(domainin.c_str());
	datain.open(filename1.c_str());
	datain >> Kn;				datain.ignore(200,'\n');
	datain >> Kt;				datain.ignore(200,'\n');
	datain >> Gn;				datain.ignore(200,'\n');
	datain >> Gt;				datain.ignore(200,'\n');
	datain >> Gv;				datain.ignore(200,'\n');
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
	datain >> starttag;  		datain.ignore(200,'\n');
	datain >> numbershapes;		datain.ignore(200,'\n');
	datain >> numberintervals; 	datain.ignore(200,'\n');
	datain.close();
	cout << "read data \n";

    Dict P;
    for (int i=0;i<numberintervals;i++)
    {
    	P.Set(/*Tag*/i*numbershapes,"Kn Kt Gn Gt Gv Mu",Kn,Kt,Gn,Gt,Gv, Mu);
        for (int j=1; j<numbershapes; j++)
        	{
        		P.Set(/*Tag*/-i*numbershapes-j,"Kn Kt Gn Gt Gv Mu",Kn,Kt,Gn,Gt, Gv, Mu);
        	}
    }
    particles.SetProps(P);
    cout << "set parameters \n";
    Vec3_t g(0.0,0.0,-9.8);
    particles.Alpha = Alpha;
    for (size_t i=0;i<particles.Particles.Size();i++)
        {
            particles.Particles[i]->Ff = particles.Particles[i]->Props.m*g;
            if (particles.Particles[i]->Tag <-999)
            	{
            		particles.Particles[i]->FixVeloc();
            	}
        }

    particles.Solve  (/*tf*/tf, /*dt*/dt, /*dtOut*/dtout, NULL, NULL, /*filekey*/filekey.c_str(),/*Visit visualization*/visualization,/*N_proc*/processornumber, /*kinematic energy*/Kinematicenergy);
    particles.Save (domainout.c_str());
    cout << "saved domain \n";
    return 0;
}
MECHSYS_CATCH
