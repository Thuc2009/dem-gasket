#ifndef MECHSYS_DEM_DOMAIN_H
#define MECHSYS_DEM_DOMAIN_H

// Std lib
#include <cmath>
#include <stdlib.h> // for M_PI
#include <iostream>
#include <fstream>
#include <set>
#include <list>
#include <utility>


// Voro++
//#include "src/voro++.cc"
#include "voro++.hh"

// MechSys
#include <mechsys/dem/interacton.h>
#include <mechsys/util/array.h>
#include <mechsys/util/util.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/mesh/mesh.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/stopwatch.h>
#include <mechsys/util/tree.h>

namespace DEM
{
struct MtData;
class Domain
{
public:
	//typedefs
	typedef void (*ptFunt_t) (Domain & Dom, void * UserData);
	// Constructor
	Domain (void * UserData =NULL)
	// Destructor
	~Domain();
	// Particle generation
	void AddDodeca(int Tag, Vec3_t const X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);
	void AddIcosa(int Tag, Vec3_t const X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);
	void AddOcta(Vec3_t const X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);
};

inline void AddOcta(int Tag, Vec3_t const X, double R, double L, double rho, double Angle, Vec3_t * Axis)
{
//Vertices
	Array<Vec3_t> V(6);
	V[0]=-L/2, -L/2, 0.0;
	V[1]=L/2, -L/2,0.0;
	V[2]=L/2, L/2, 0.0;
	V[3]=-L/2, L/2, 0.0;
	V[4]=0.0, 0.0, L/sqrt(2.0);
	V[5]=0.0, 0.0, -L/sqrt(2.0);
//Edges
	Array<Array<int> > E(12);
	for (size_t i=0;i<12;i++) E[i].Resize(2);
	E[0]=0,1;
	E[1]=1,2;
	E[2]=2,3;
	E[3]=3,0;
	E[4]=0,4;
	E[5]=1,4;
	E[6]=2,4;
	E[7]=3,4;
	E[8]=0,5;
	E[9]=1,5;
	E[10]=2,5;
	E[11]=3,5;
//Faces
	Array<Array<int> >F(8);
	for (size_t i=0;i<8;i++) F[i].resize[3];
	F[0]=0,1,4;
	F[1]=1,2,4;
	F[2]=2,3,4;
	F[3]=3,0,4;
	F[4]=0,1,5;
	F[5]=1,2,5;
	F[6]=2,3,5;
	F[7]=3,0,5;
//Rotation
	bool ThereisanAxis =true;
	if (Axis==NULL)
	{
        Angle   = (1.0*rand())/RAND_MAX*2*M_PI;
        Axis = new Vec3_t((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
        ThereisanAxis = false;
	}
    Quaternion_t q;
    NormalizeRotation (Angle,(*Axis),q);
    for (size_t i=0; i<V.Size(); i++)
    {
        Vec3_t t;
        Rotation (V[i],q,t);
        V[i] = t+X;
    }
    // add particle
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));

    // clean up
    if (!ThereisanAxis) delete Axis;
    q(0) = 1.0;
    q(1) = 0.0;
    q(2) = 0.0;
    q(3) = 0.0;
    q = q/norm(q);
    Particles[Particles.Size()-1]->Q          = q;
    Particles[Particles.Size()-1]->Props.V    = sqrt(2.0)*L*L*L/12.0;
    Particles[Particles.Size()-1]->Props.m    = rho*sqrt(2.0)*L*L*L/12.0;
    Particles[Particles.Size()-1]->I          = L*L, L*L, L*L;
    Particles[Particles.Size()-1]->I         *= Particles[Particles.Size()-1]->Props.m/20.0;
    Particles[Particles.Size()-1]->x          = X;
    Particles[Particles.Size()-1]->Ekin       = 0.0;
    Particles[Particles.Size()-1]->Erot       = 0.0;
    Particles[Particles.Size()-1]->Dmax       = sqrt(3.0*L*L/8.0)+R;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index      = Particles.Size()-1;
}
}
