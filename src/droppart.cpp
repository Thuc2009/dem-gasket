#include <mechsys/dem/domain.h>
#include <stdlib.h>
using namespace std;
using namespace DEM;
void Setup (DEM::Domain & dom, void * UD)
{
	std::cout <<"all the time: hi" << dom.Time << "  dt=" << std::endl;
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
	int numberparts;
	int startpart;
	string filekey;
	int visualization;
	int processornumber;
	double roundratio;
	DEM::Domain particles;
	//DEM::Domain particlesbefore;
	Array<Vec3_t> posbefore;
	Vec3_t xmin;
	Vec3_t xmax;
	Vec3_t xminarrange;
	Vec3_t xmaxarrange;
	Vec3_t xminup;
	Vec3_t xmaxup;
	Vec3_t xmindown;
	Vec3_t xmaxdown;
	Vec3_t move;
	Vec3_t partdepth;
	String file;
	string datafile="droppart.drp";
	string domainin;
	string domainout;
	ifstream datain;
	datain.open(datafile.c_str());
	datain >> domainin; 				datain.ignore(200,'\n');
	datain >> domainout;				datain.ignore(200,'\n');
	datain.close();
	particles.Load(domainin.c_str());
	datain.open("droppart.par");
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
	datain >> numberparts;		datain.ignore(200,'\n');
	datain >> startpart;		datain.ignore(200,'\n');
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
    Vec3_t g(0.0,0.0,-9.8);
    particles.Alpha = Alpha;
    for (size_t i=0;i<particles.Particles.Size();i++)
        {
            particles.Particles[i]->Ff = particles.Particles[i]->Props.m*g;
            posbefore.Push(particles.Particles[i]->x);
        }
    cout << "set parameters \n";
	particles.BoundingBox(xmin, xmax);
	Vec3_t z(0,0,1);
	partdepth = dot(xmax-xmin,z)/numberparts*z;
	cout << partdepth<< "\n";
	xminarrange = xmin;
	xmaxarrange=xmax;
	DEM::Domain partparticles[numberparts];
    if (numberparts>0)
    	{
    		for (int i =startpart; i<numberparts;i++)
    			{
    				xminarrange=xmin+i*partdepth;
    				xmaxarrange=xmin+(i+1)*partdepth;
    				cout << xminarrange << "\n";
    				cout << xmaxarrange << "\n";
    				for (size_t j=0; j<particles.Particles.Size(); j++)
    					{
    						if ((particles.Particles[j]->x(2)>xminarrange(2))and(particles.Particles[j]->x(2) <xmaxarrange(2)))
    							{
    								partparticles[i].Particles.Push(particles.Particles[j]);
    							}
    					}
    				if (i>0)
    					{
        					partparticles[i-1].BoundingBox(xmindown, xmaxdown);
        					partparticles[i].BoundingBox(xminup,xmaxup);
        					move = dot(xmaxdown-xminup,z)*z;
        					for (size_t j=0; j<partparticles[i].Particles.Size(); j++)
        						{
        							partparticles[i].Particles[j]->Translate(move);
        						}
    					}
    				partparticles[i].GenBoundingBox(starttag,0.01,1.2,false);
    				for (int j=starttag;j>starttag-6;j--)
    				    {
    				        partparticles[i].GetParticle(j)->FixVeloc();
    				    }
    				file.Printf("%s_%04d", filekey.c_str(), i);
                    cout << "Part: "<< i << "particles: "<< partparticles[i].Particles.Size()<< "\n";
                    if (partparticles[i].Particles.Size()>0)
                    	{
                    		partparticles[i].Solve(/*tf*/ tf/numberparts, /*dt*/dt, /*dtOut*/dtout, NULL, NULL, /*filekey*/file.c_str(),/*Visit visualization*/visualization,/*N_proc*/processornumber, /*kinematic energy*/Kinematicenergy);
                    		particles.Save(file.c_str());
                    		particles.WriteXDMF(file.c_str());// export to draw visual results
                    	}
    			}
    	}
    particles.Solve  (/*tf*/tf, /*dt*/dt, /*dtOut*/dtout, NULL, NULL, /*filekey*/filekey.c_str(),/*Visit visualization*/visualization,/*N_proc*/processornumber, /*kinematic energy*/Kinematicenergy);
    particles.Save (domainout.c_str());
    cout << "solve domain \n";

    //for (size_t i=0;i<particles.Particles.Size();i++)
    //    {
    //       //if (norm(particles.Particles[i]->x-particlesbefore.Particles[i]->x)>0.1)
    //	    if (norm(particles.Particles[i]->x-posbefore[i])>0.1)
    //        {
    //        	count1++;
    //       }
    //        if (norm(particles.Particles[i]->x-posbefore[i])> 2*particles.Particles[i]->Props.R)
    //        {
    //        	count2++;
    //        }
    //    }
    //cout << count1 << " " << count2;
    return 0;
}
MECHSYS_CATCH
