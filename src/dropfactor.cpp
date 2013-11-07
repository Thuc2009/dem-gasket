#include <mechsys/dem/domain.h>
using namespace std;
using namespace DEM;
void Report (DEM::Domain & particles, void *UD)
{
	//UserData & dat = (*static_cast<UserData *>(UD));
	String fn;
	fn.Printf    ("%s_%s_%04d", particles.FileKey.CStr(), "bf",particles.idx_out);
//	particles.WriteVTKContacts(fn.CStr());
//	particles.WriteBF(fn.CStr());
	particles.Save (fn.c_str());
	//dat.pressure = 0;
}
void AddFunnel(DEM::Domain & dom, int Tag,  Vec3_t & X0, Vec3_t & X1, double Lx0, double Ly0, double Lx1, double Ly1, double rho, double R);
int main(int argc, char **argv) try
{
	string filename;
	string filename1;
	if (argc<2)
		{
			filename = "dropfactor.drp";
			filename1 = "dropfactor.par";
		}
	else if (argc <3)
		{
			filename = argv[1];
			filename1 = "dropfactor.par";
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
	//DEM::Domain particlesbefore;
	//Array<Vec3_t> posbefore;
	Vec3_t Xmin;
	Vec3_t Xmax;
	Array <int> delparbefore;
	Array <int> delparafter;
	double size[3];
	string domainin;
	string domainout;
	string boundary;
	ifstream datain;
	datain.open(filename.c_str());
	datain >> domainin; 				datain.ignore(200,'\n');
	datain >> domainout;				datain.ignore(200,'\n');
	datain >> boundary;					datain.ignore(200,'\n');
	if (boundary =="cylinder")
		{
			for (int i=0; i<2; i++)
				{
					datain>>size[i];
				}
		}
	else if (boundary =="box")
		{
			for (int i=0; i<3; i++)
				{
					datain>>size[i];
				}
		}
	datain.close();
	particles.Load(domainin.c_str());
	particles.BoundingBox(Xmin,Xmax);
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
	Vec3_t X3 = (Xmax+Xmin)/2+OrthoSys::e2*dot((Xmax-Xmin)/2,OrthoSys::e2);
	Vec3_t X2 = (Xmax+Xmin)/2-OrthoSys::e2*dot((Xmax-Xmin)/2,OrthoSys::e2);
	Vec3_t X1 = X2 -OrthoSys::e2*dot((Xmax-Xmin)/2,OrthoSys::e2);
	Vec3_t X0 = X1 -Vec3_t(0.,0.,size[2]);
	if (boundary == "cylinder")
		{
			// generate bounding cylinder
			double Lx1 = dot(Xmax-Xmin,Vec3_t(OrthoSys::e0));
			double Ly1 = dot(Xmax-Xmin,Vec3_t(OrthoSys::e1));
			double r= pow(Lx1*Lx1+Ly1*Ly1,0.5);
			particles.AddCylinder(starttag,X2,r,X3,r,0.05,1.);
			particles.AddPlane(starttag-1,X3,0.05,2.4*r,2.4*r,1.);
			// generate new cylinder
			particles.AddCylinder(2*starttag,X0,size[0],X1,size[0],0.05,1.);
			particles.AddPlane(2*starttag-1,X0,0.05,1.2*size[0],1.2*size[0],1.);
			// generate new cylinder funnel
			particles.AddCylinder(3*starttag,X1,size[0],X2,r,0.05,1.);
			particles.GetParticle(starttag)->FixVeloc();
			particles.GetParticle(starttag-1)->FixVeloc();
			particles.GetParticle(2*starttag)->FixVeloc();
			particles.GetParticle(2*starttag-1)->FixVeloc();
			particles.GetParticle(3*starttag)->FixVeloc();
		}
	else if (boundary=="box")
		{
			// generate bounding box
			particles.GenBoundingBox(starttag,0.01,1.2,false);
			for (int i=starttag;i>starttag-5;i--)
		    	{
		        	particles.GetParticle(i)->FixVeloc();
		    	}
			delparbefore.Push(starttag-5);
			particles.DelParticles(delparbefore);
			//generate new box
			particles.GenOpenBox(2*starttag,size[0],size[1],size[2],0.05,1.0);
			Vec3_t X= (X0+X1)/2-(Xmax+Xmin)/2;						// for vertical funnel
//			Vec3_t X=-2*(Xmax-Xmin);						// for inclined funnel
			for (int i =0; i<5;i++)
				{
					particles.GetParticle(2*starttag-i)->Translate(X);
					particles.GetParticle(2*starttag-i)->FixVeloc();
				}
			//generate funnel connecting two boxes

			double Lx1 = dot(Xmax-Xmin,Vec3_t(OrthoSys::e0));
			double Ly1 = dot(Xmax-Xmin,Vec3_t(OrthoSys::e1));
			AddFunnel(particles,3*starttag,X0,X1, size[0], size[1], Lx1 , Ly1, 1.0, 0.05);
			particles.GetParticle(3*starttag)->FixVeloc();

		}
//	particles.Save (domainout.c_str());
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
            //posbefore.Push(particles.Particles[i]->x);
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
//    particles.Solve  (/*tf*/tf, /*dt*/dt, /*dtOut*/dtout, NULL, &Report, /*filekey*/filekey.c_str(),/*Visit visualization*/visualization,/*N_proc*/processornumber, /*kinematic energy*/Kinematicenergy);
    particles.Solve  (/*tf*/tf, /*dt*/dt, /*dtOut*/dtout, NULL, &Report, /*filekey*/filekey.c_str(),/*Visit visualization*/visualization,/*N_proc*/processornumber, /*kinematic energy*/Kinematicenergy);
    for (size_t i; i<particles.Particles.Size();i++)
    	{
    		if (particles.Particles[i]->Tag <-999)
    			{
    				delparafter.Push(particles.Particles[i]->Tag);
    			}
    	}
    particles.DelParticles(delparafter);
    particles.Save (domainout.c_str());
    cout << "saved domain \n";

//    for (size_t i=0;i<particles.Particles.Size();i++)
//        {
//            //if (norm(particles.Particles[i]->x-particlesbefore.Particles[i]->x)>0.1)
//    	    if (norm(particles.Particles[i]->x-posbefore[i])>0.1)
//            {
//            	count1++;
//            }
//            if (norm(particles.Particles[i]->x-posbefore[i])> 2*particles.Particles[i]->Props.R)
//            {
//            	count2++;
//            }
//        }
//    cout << count1 << " " << count2;

    return 0;
}
MECHSYS_CATCH

inline void AddFunnel(DEM::Domain & dom, int Tag,  Vec3_t & X0, Vec3_t & X1, double Lx0, double Ly0, double Lx1, double Ly1, double rho, double R)
{
     Array<Vec3_t>        V(8);
     Array<Array <int> >  E(12);
     Array<Array <int> >  F(4);

     V[0]=X0-Lx0/2*OrthoSys::e0-Ly0/2*OrthoSys::e1;
     V[1]=X0+Lx0/2*OrthoSys::e0-Ly0/2*OrthoSys::e1;
     V[2]=X0+Lx0/2*OrthoSys::e0+Ly0/2*OrthoSys::e1;
     V[3]=X0-Lx0/2*OrthoSys::e0+Ly0/2*OrthoSys::e1;
     E[0].Push(0); E[0].Push(1);
     E[1].Push(1); E[1].Push(2);
     E[2].Push(2); E[2].Push(3);
     E[3].Push(3); E[3].Push(0);

     V[4]=X1-Lx1/2*OrthoSys::e0-Ly1/2*OrthoSys::e1;
     V[5]=X1+Lx1/2*OrthoSys::e0-Ly1/2*OrthoSys::e1;
     V[6]=X1+Lx1/2*OrthoSys::e0+Ly1/2*OrthoSys::e1;
     V[7]=X1-Lx1/2*OrthoSys::e0+Ly1/2*OrthoSys::e1;
     E[4].Push(4); E[4].Push(5);
     E[5].Push(5); E[5].Push(6);
     E[6].Push(6); E[6].Push(7);
     E[7].Push(7); E[7].Push(4);

     E[8].Push(0); E[8].Push(4);
     E[9].Push(1); E[9].Push(5);
     E[10].Push(2); E[10].Push(6);
     E[11].Push(3); E[11].Push(7);
     F[0].Push(0);F[0].Push(1);F[0].Push(5);F[0].Push(4);
     F[1].Push(1);F[1].Push(2);F[1].Push(6);F[1].Push(5);
     F[2].Push(2);F[2].Push(3);F[2].Push(7);F[2].Push(6);
     F[3].Push(3);F[3].Push(0);F[3].Push(4);F[3].Push(7);
     //std::cout << "2" << std::endl;
     dom.Particles.Push(new DEM::Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
     //std::cout << "3" << std::endl;
     dom.Particles[dom.Particles.Size()-1]->Q          = 1.0,0.0,0.0,0.0;
     dom.Particles[dom.Particles.Size()-1]->Props.V    = (2*(Lx0+Lx1)*pow(pow(dot(X1-X0,OrthoSys::e2),2.0)+pow((Ly1-Ly0)/2.,2.0),0.5)
    													 +2*(Ly0+Ly1)*pow(pow(dot(X1-X0,OrthoSys::e2),2.0)+pow((Lx1-Lx0)/2.,2.0),0.5))*R;
     dom.Particles[dom.Particles.Size()-1]->Props.m    = rho*dom.Particles[dom.Particles.Size()-1]->Props.V;
     dom.Particles[dom.Particles.Size()-1]->I          = 1.0,1.0,1.0;
     dom.Particles[dom.Particles.Size()-1]->I         *= dom.Particles[dom.Particles.Size()-1]->Props.m;
     dom.Particles[dom.Particles.Size()-1]->x          = (X0+X1)/2;
     dom.Particles[dom.Particles.Size()-1]->Ekin       = 0.0;
     dom.Particles[dom.Particles.Size()-1]->Erot       = 0.0;
     dom.Particles[dom.Particles.Size()-1]->Dmax       = sqrt(Lx1*Lx1+Ly1*Ly1)+R;
     dom.Particles[dom.Particles.Size()-1]->PropsReady = true;
     dom.Particles[dom.Particles.Size()-1]->Index      = dom.Particles.Size()-1;

     dom.Particles[dom.Particles.Size()-1]->Position((X0+X1)/2);
}
