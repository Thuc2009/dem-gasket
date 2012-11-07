#include <mechsys/dem/domain.h>

void Setup (DEM::Domain & dom, void * UD)
{
	std::cout <<"all the time: hi" << dom.Time << "  dt=" << std::endl;
}

void Report (DEM::Domain & dom, void * UD)
{
	dom.Save("mysim");
	std::cout <<"     only after time out: hi" << dom.Time << "  dt=" << std::endl;
}

int main(int argc, char **argv) try
{
    // domain
    DEM::Domain dom;
    for (int k =0; k<20 ;k++)
    {
        for (int i=0; i <4; i++)
        	{
        		for (int j=0; j<4; j++)
        		{
            		if ((i+j+k) % 3 ==0)
            				{
            					dom.AddCube(-i-j,Vec3_t(i, j, k),0.05 /*R*/,0.5*(rand()/RAND_MAX+1.0)/*Length*/,3.2 /*rho*/,M_PI/3.0,&OrthoSys::e0 );
            				}
            			else if ((i+j+k)%3==1)
            				{
            					dom.AddTetra(-i-j,Vec3_t(i, j, k),0.05 /*R*/,0.5*(rand()/RAND_MAX+1.0)/*Length*/,3.2 /*rho*/,M_PI/3.0,&OrthoSys::e0 );
            				}
            			else
            				{
            					dom.AddSphere(i+j,Vec3_t(i, j, k),0.25*(rand()/RAND_MAX+1.0),3.2 /*rho*/);
            				}
        		}
        	}
    }

    dom.GenBoundingBox(-100, 0.05, 1.0,false);
    dom.GetParticle(-100)->FixVeloc();
    dom.GetParticle(-101)->FixVeloc();
    dom.GetParticle(-102)->FixVeloc();
    dom.GetParticle(-103)->FixVeloc();
    dom.GetParticle(-104)->FixVeloc();
    dom.GetParticle(-105)->FixVeloc();

    //Add the gravity
    Vec3_t g(0.0,0.0,-9.8);
    for (size_t i=0;i<dom.Particles.Size();i++)
    {
        dom.Particles[i]->Ff = dom.Particles[i]->Props.m*g;
        dom.Particles[i]->Props.Kn = 5.0e4;    // Value for the normal stiffness
        dom.Particles[i]->Props.Kt = 5.0e4;    // Value for the tangential stiffness.
        dom.Particles[i]->Props.Mu = 0.2;    // Value for friction coefficient
    }

    dom.Solve     (/*tf*/10.0, /*dt*/1.0e-4, /*dtOut*/0.1, &Setup, &Report, /*filekey*/"dem",/*Visit visualization*/2,/*N_proc*/3);
    dom.Save("dem_sample");
    return 0;
}
MECHSYS_CATCH
