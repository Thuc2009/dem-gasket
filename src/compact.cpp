#include <mechsys/dem/domain.h>
using namespace std;
using namespace DEM;

struct UserData
{
	double pressure;
};

void Report (DEM::Domain & particles, void *UD)
{
	//UserData & dat = (*static_cast<UserData *>(UD));
	String fn;
	fn.Printf    ("%s_%s_%04d", particles.FileKey.CStr(), "bf",particles.idx_out);
	particles.WriteVTKContacts(fn.CStr());
	particles.WriteBF(fn.CStr());
	//dat.pressure = 0;
}

int main(int argc, char **argv) try
{
	string filename;
	string filename1;
	if (argc<2)
		{
			filename = "compact.drp";
			filename1 = "compact.par";
		}
	else if (argc <3)
		{
			filename = argv[1];
			filename1 = "compact.par";
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
	double Mu;
	double dt;
	double tf;
	double dtout;
	double Alpha;
	double Kinematicenergy;
	double starttag;
	int numbershapes;
	int numberintervals;
	string filekey;
	int visualization;
	int processornumber;
	double roundratio;
	//UserData dat;
	//DEM::Domain particles(&dat);
	DEM::Domain particles;
	//DEM::Domain particlesbefore;
	Array<Vec3_t> posbefore;
	string domainin;
	string domainout;
	int tag;
	String axis;
	double force;
	ifstream datain;
	datain.open(filename.c_str());
	datain >> domainin; 				datain.ignore(200,'\n');
	datain >> domainout;				datain.ignore(200,'\n');
	datain >> tag;						datain.ignore(200,'\n');
	datain >> axis;						datain.ignore(200,'\n');
	datain >> force;					datain.ignore(200,'\n');
	datain.close();
	particles.Load(domainin.c_str());
	datain.open(filename1.c_str());
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
	datain >> starttag;			datain.ignore(200,'\n');
	datain >> numbershapes;		datain.ignore(200,'\n');
	datain >> numberintervals; 	datain.ignore(200,'\n');
	datain.close();
	cout << "data reading \n";
    Vec3_t g(0.0,0.0,-0.098);
    //for (size_t i=0;i<particles.Particles.Size();i++)
    //    {
    //      particles.Particles[i]->Ff = particles.Particles[i]->Props.m*g;
    //        posbefore.Push(particles.Particles[i]->x);
    //    }
	if (starttag < 0)
		{
			Array <int> delpar;
			double z =0.;
			int top =0;
			for (size_t i=0; i<particles.Particles.Size(); i++)
				{
					if (particles.Particles[i]->Tag <-999)
						{
							delpar.push_back(particles.Particles[i]->Tag);
						}
					else if (particles.Particles[i]->x[2]>z)	// find the top particles
						{
							z=particles.Particles[i]->x[2];
							top = i;
						}
				}
			delpar.push_back(top);
			particles.DelParticles(delpar);
			particles.GenBoundingBox(starttag,0.01,1.2,false);
			for (int i=starttag; i>starttag-6; i--)
			    {
			        particles.GetParticle(i)->FixVeloc();
			    }
		}
	if (axis=="z")
	{
		particles.GetParticle(starttag+tag)->vzf = false;
		particles.GetParticle(starttag+tag)->Ff=Vec3_t(0,0,force);
	}
    Dict P;
    for (size_t i=0;i<numberintervals;i++)
    {
        P.Set(/*Tag*/i*numbershapes,"Kn Kt Gn Gt Mu",Kn,Kt,Gn,Gt,Mu);
        for (int j=1; j<numbershapes; j++)
        	{
        		P.Set(/*Tag*/-i*numbershapes-j,"Kn Kt Gn Gt Mu",Kn,Kt,Gn,Gt,Mu);
        	}
    }
    particles.SetProps(P);
    cout << "set parameters \n";
    particles.Alpha = Alpha;
    int time=0;
    particles.Solve  (/*tf*/tf, /*dt*/dt, /*dtOut*/dtout, NULL, &Report, /*filekey*/filekey.c_str(),/*Visit visualization*/visualization,/*N_proc*/processornumber, /*kinematic energy*/Kinematicenergy);
    particles.Save (domainout.c_str());
    cout << "solve domain \n";
	return 0;
}
MECHSYS_CATCH
