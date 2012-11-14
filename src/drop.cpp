#include <mechsys/dem/domain.h>
using namespace std;
using namespace DEM;
void Report (DEM::Domain & particles, void * UD)
{
	particles.Save ("drop");
}
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
	int starttag;
	int numbershapes;
	int numberintervals;
	string filekey;
	int visualization;
	int processornumber;
	double roundratio;
	DEM::Domain particles;
	//DEM::Domain particlesbefore;
	Array<Vec3_t> posbefore;
	string datafile="drop.drp";
	string domainin;
	string domainout;
	ifstream datain;
	datain.open(datafile.c_str());
	datain >> domainin; 				datain.ignore(200,'\n');
	datain >> domainout;				datain.ignore(200,'\n');
	datain.close();
	particles.Load(domainin.c_str());
	datain.open("drop.par");
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
	datain >> starttag;  		datain.ignore(200,'\n');
	datain >> numbershapes;		datain.ignore(200,'\n');
	datain >> numberintervals; 	datain.ignore(200,'\n');
	datain.close();
	//particlesbefore.Load(domainin.c_str());
	cout << "read data \n";
	if (starttag<0)
		{
			particles.GenBoundingBox(starttag,0.01,1.2,false);
			for (int i=starttag;i>starttag-6;i--)
			    {
			        particles.GetParticle(i)->FixVeloc();
			    }
		}
    Dict P;
    for (int i=0;i<numberintervals;i++)
    {
        P.Set(/*Tag*/i*numbershapes,"Kn Kt Gn Gt Mu",Kn,Kt,Gn,Gt,Mu);
        for (int j=1; j<numbershapes; j++)
        	{
        		P.Set(/*Tag*/-i*numbershapes-j,"Kn Kt Gn Gt Mu",Kn,Kt,Gn,Gt,Mu);
        	}
    }
    particles.SetProps(P);
    cout << "set parameters \n";
    Vec3_t g(0.0,0.0,-9.8);
    particles.Alpha = Alpha;
    for (size_t i=0;i<particles.Particles.Size();i++)
        {
            particles.Particles[i]->Ff = particles.Particles[i]->Props.m*g;
            posbefore.Push(particles.Particles[i]->x);
        }

    //for (size_t i=0;i<particles.Particles.Size()-1;i++)
    //        {
    //			for (size_t j=i+1;j<particles.Particles.Size();j++)
    //	        	{
    //					if (norm(particles.Particles[i]->x-particles.Particles[j]->x)<(particles.Particles[i]->Props.R+particles.Particles[j]->Props.R-0.001))
    //					{
    //						count1++;
    //						cout << i << " " <<particles.Particles[i]->Tag << " : " << j << " " << particles.Particles[j]->Tag << "\n";
    // 					}
    //	        	}
    //       }
    //cout << count1 <<"\n";
    //count1=0;
    particles.Solve  (/*tf*/tf, /*dt*/dt, /*dtOut*/dtout, NULL, &Report, /*filekey*/filekey.c_str(),/*Visit visualization*/visualization,/*N_proc*/processornumber, /*kinematic energy*/Kinematicenergy);
    particles.Save (domainout.c_str());
    cout << "solve domain \n";

    for (size_t i=0;i<particles.Particles.Size();i++)
        {
            //if (norm(particles.Particles[i]->x-particlesbefore.Particles[i]->x)>0.1)
    	    if (norm(particles.Particles[i]->x-posbefore[i])>0.1)
            {
            	count1++;
            }
            if (norm(particles.Particles[i]->x-posbefore[i])> 2*particles.Particles[i]->Props.R)
            {
            	count2++;
            }
        }
    cout << count1 << " " << count2;
    return 0;
}
MECHSYS_CATCH
