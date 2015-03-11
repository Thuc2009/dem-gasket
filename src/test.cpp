// Std lib
#include <cmath>
#include <stdlib.h> // for M_PI
#include <iostream>
#include <fstream>
#include <set>
#include <list>
#include <utility>
#include <boost/date_time/posix_time/posix_time.hpp> // for time to milliseconds

// Voro++
#include "voro++.hh"

// MechSys

#include <mechsys/dem/domain.h>
#include <mechsys/dem/particle.h>
#include "./gsd.h"
// using namespace
using namespace std;
using namespace DEM;
// this namespace
int main(int argc, char **argv)
{
	DEM::Domain dom;
	DEM::Particle * Pa = new DEM::Particle;
	double R=0.01;
	double rho=1;
	Pa->ConstructFromJson(-1,"octahedron.msh",R,rho,1.0);
	dom.Particles.Push(Pa);
	Pa->ConstructFromJson(-1,"dodecahedron.msh",R,rho,1.0);
	dom.Particles.Push(Pa);
	Pa->ConstructFromJson(-1,"icosahedron.msh",R,rho,1.0);
	dom.Particles.Push(Pa);
	for (int i=0; i<3; i++)
	{
		dom.Particles[i]->Shrink(1.0/dom.Particles[i]->Dmax);
		dom.Particles[i]->x= Vec3_t (i,0,0);
	}
	dom.WriteXDMF("test");
	dom.Save("test.dom");
}
