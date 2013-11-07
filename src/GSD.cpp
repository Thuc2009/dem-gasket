// Std lib
#include <cmath>
#include <stdlib.h> // for M_PI
#include <iostream>
#include <fstream>
#include <set>
#include <list>
#include <utility>

// Voro++
#include "voro++.hh"

// MechSys

#include <mechsys/dem/domain.h>
// using namespace
using namespace std;
// this namespace
namespace DEM
{
struct SequentialFace
{
	int points[3];
	int point;
	int faceuse;
	double constrictionsize;
	Vec3_t normal;
	Vec3_t constrictioncentre;
};
struct SequentialTetrahedron
{
	int points[4];
	int faces[4];
	double poreradius;
	Vec3_t porecentre;
};
struct File
{
	string filein;
	string fileout;
	string parameterfile;
};
struct Interval
{
	bool ability;
	int firstparticle;														// particle number
	int lastparticle;
	int usingparticle;
	double begin;															// fraction
	double begindiameter;
	double diameter;														//mean diameter;
	double end;
	double enddiameter;
	double fraction;
	double mass;
	double numberparticles;
	double specialgravity;
	double volumeparticles;
	double shaperatios[10]; 	// 0-sphere 1- cube 2- tetrahedron 3-rectangular box 9-user defined
	int shapestates[10][4];	// 0 - numberparticles, 1- begin ID, 2 - end ID, 3 - current using
};
struct Gsd
{
	string Soilname;
	string arrangement; 	// particle arrangement 0 - customized, 1- layer, 2 - discrete
	string boundaryfile;
	double density;
	double mass;
	double numberparticles;
	double porosity;
	double specialgravity;	// special gravity 0 - provided by intervals
	double volume;
	int maxfraction; 		// max percentage per intervals
	int numberintervals;
	int numberinters;
	int numbershapes;
	double shapemass[10];	//ratio between mass of shape and basic size
	double rectangularboxratios[3];
	vector <Interval> intervals;
	vector <Interval> inters;
};
struct SequentialPackingData
{
	vector <SequentialFace> faces;
	vector <SequentialTetrahedron> tetrahedrons;
	vector <bool> particleuses;
	vector <int> boundary;
	vector <int> threadnumbers;
	Vec3_t localroot;
	Vec3_t localsystem[3];
	bool check;
	bool checkboundary;
	bool checkoverlap;
	bool checkradius;
	bool facereuse;
	bool particlereuse;
	double approximation;
	double boundaryfactors[10];
	double boundarysizes[10];
	int facereusenumber;
	int numberopenfaces;
	int numberfaces;
	int numberprocessors;
	int numberunusedparticles;
	int overlappingpoint;
	int usingface;
	int usinginter;
	int usingparticle;
	int temproraryparticle;
	string boundarytype;
	DEM::Gsd gsd;
	DEM::Domain specimen;
};

// main classes
class Soil
{
public:
    // typedefs
    typedef void (*ptFun_t) (Soil & soil, void * UserData);

    // Constructor
    Soil(void * UserData=NULL);

    // Destructor
    ~Soil();
    // List of functions

    double ShapeMass(DEM::SequentialPackingData&packinfo, int shape, double size, double density);
    void AddParticle( DEM::SequentialPackingData&packinfo, int inter, int shape, double size, double density, Vec3_t position, double roundness=0.005);
    void BasicTetrahedron (DEM::SequentialPackingData&packinfo, int p[4]);
    void CheckBoundary(DEM::SequentialPackingData&packinfo, int p3);
    void CheckOverlap( DEM::SequentialPackingData&packinfo, int p3);
    void CloseSequentialFace(DEM::SequentialPackingData&packinfo, int face);
    void CompleteGsdData(DEM::SequentialPackingData&packinfo);
    void CreateSequentialFace(DEM::SequentialPackingData&packinfo, int p[3],int p3, bool sort=true);
    void CreateSequentialTetrahedron(DEM::SequentialPackingData&packinfo, int face, int p3);
    void DrawBoundary(DEM::SequentialPackingData&packinfo);
    void EstablishLocalSystem(DEM::SequentialPackingData&packinfo, int face);
    void FrozenTag(DEM::SequentialPackingData&packinfo, int tag);
    void PrepareList(DEM::SequentialPackingData&packinfo, bool randomness=false);
	void PutSequentialParticle(DEM::SequentialPackingData&packinfo, int face, int p3);
    void ReadGsd(DEM::SequentialPackingData&packinfo,string FileKey="sand");
    void SequentialPacking(DEM::SequentialPackingData&packinfo);
    void TrySequentialParticle(DEM::SequentialPackingData&packinfo, int face, int p3);
    void UseParticle(DEM::SequentialPackingData&packinfo, int p3);
};

// Constructor & Destructor
inline Soil::Soil(void * UD)
{}

inline Soil::~Soil()
{}



inline double Soil::ShapeMass(DEM::SequentialPackingData&packinfo, int shape, double size, double density)
{
	double mass;
	switch (shape)
	{
	case 0:
		mass = density * acos(0.) * pow(size,3.) / 3.;
		break;
	case 1:
		mass = density * pow(size,3.);
		break;
	case 2:
		mass = density * pow(size,3.) / 3.;
		break;
	case 3:
		mass = density * pow(pow(2.,0.5)*size/(1+packinfo.gsd.rectangularboxratios[1]),3.) *packinfo.gsd.rectangularboxratios[1] * packinfo.gsd.rectangularboxratios[2];
		break;
	default :
		break;
	}
	return (mass);
}
// Functions
inline void Soil::AddParticle( DEM::SequentialPackingData&packinfo, int inter, int shape, double size, double density, Vec3_t position, double roundness)
{
	switch (shape)
	{
	case 0:
		packinfo.specimen.AddSphere(packinfo.gsd.numbershapes*inter,position,size/2,density);
		break;
	case 1:
		packinfo.specimen.AddCube(-packinfo.gsd.numbershapes*inter-shape,position,roundness,size,density);
		packinfo.specimen.Particles[packinfo.specimen.Particles.Size()-1]->Erode(roundness);
		break;
	case 2:
		packinfo.specimen.AddTetra(-packinfo.gsd.numbershapes*inter-shape,position,roundness,pow(2.,0.5)*size,density);
		packinfo.specimen.Particles[packinfo.specimen.Particles.Size()-1]->Erode(roundness);
		break;
	case 3:
		double smallsize;
		smallsize= size/(1+packinfo.gsd.rectangularboxratios[1]);
		packinfo.specimen.AddRecBox(-packinfo.gsd.numbershapes*inter-shape,position,Vec3_t (smallsize,smallsize*packinfo.gsd.rectangularboxratios[1],smallsize*gsd.rectangularboxratios[2]),roundness,density);
		packinfo.specimen.Particles[packinfo.specimen.Particles.Size()-1]->Erode(roundness);
		break;
	default:
		break;
	}
}

inline void Soil::BasicTetrahedron(DEM::SequentialPackingData&packinfo, int p[4])
{
	packinfo.numberunusedparticles = packinfo.specimen.Particles.Size();
	packinfo.numberopenfaces=0;
	cout << "Creating basic tetrahedron \n";
	double r[4];
	for (int i=0;i<4;i++)
	{
		r[i]=packinfo.specimen.Particles[p[i]]->Dmax;
	}
	packinfo.specimen.Particles[p[0]]->Position(Vec3_t(0.,0.,0.));
	packinfo.specimen.Particles[p[1]]->Position(Vec3_t(r[0]+r[1],0.,0.));
	double x = (pow(r[0]+r[1],2)+pow(r[0]+r[2],2)-pow(r[1]+r[2],2))/2/(r[0]+r[1]);
	packinfo.specimen.Particles[p[2]]->Position(Vec3_t(x, pow(pow(r[0]+r[2],2)-pow(x,2),0.5),0.));
	int f[3]={p[0],p[1],p[2]};
	CreateSequentialFace(packinfo,f,p[3],true);
	PutSequentialParticle(packinfo,0,packinfo.usingparticle-3);
	CreateSequentialTetrahedron(packinfo, 0,packinfo.usingparticle-3);
	packinfo.faces[0].faceuse=1;																											// mark used side
	packinfo.usingface=0;
	packinfo.numberfaces =4;				// or 1 ?
	packinfo.numberopenfaces =4;			// or 1 ?
	for (int i=0; i<4; i++)																							//mark used particles
		{
			UseParticle(packinfo, p[i]);
		}
}

inline void Soil::CheckBoundary( DEM::SequentialPackingData&packinfo, int p3)
{
	packinfo.checkboundary=true;
	if (packinfo.boundarytype =="File")
	{
		for (int i=0; i<packinfo.boundary.size();i++)
		{

		}
	}
	else if (packinfo.boundarytype=="Sphere")
	{
		if (norm(packinfo.specimen.Particles[p3]->x)+ packinfo.specimen.Particles[p3]->Dmax >packinfo.boundaryfactors[0])
			{
				packinfo.checkboundary=false;
			}
	}
	else if (packinfo.boundarytype=="Cube")
	{
		for (int i = 0; i < 3; i++)
			{
				if (abs(packinfo.specimen.Particles[p3]->x(i))+packinfo.specimen.Particles[p3]->Dmax > packinfo.boundaryfactors[i])
					{
						packinfo.checkboundary= false;
					}
			}
	}
	else if (packinfo.boundarytype== "Cylinder")
	{
		if (abs(packinfo.specimen.Particles[p3]->x(2))+packinfo.specimen.Particles[p3]->Dmax > packinfo.boundaryfactors[1])
			{
				packinfo.checkboundary=false;
			}
		else if (abs(pow(pow(packinfo.specimen.Particles[p3]->x(0),2.)+pow(packinfo.specimen.Particles[p3]->x(1),2.),0.5))+packinfo.specimen.Particles[p3]->Dmax > packinfo.boundaryfactors[0])
			{
				packinfo.checkboundary=false;
			}
	}
	else if (packinfo.boundarytype=="Box")
	{
		for (int i = 0; i < 3; i++)
			{
				if (abs(packinfo.specimen.Particles[p3]->x(i))+packinfo.specimen.Particles[p3]->Dmax > packinfo.boundaryfactors[i])
					{
						packinfo.checkboundary= false;
					}
			}
	}
}

inline void Soil::CheckOverlap( DEM::SequentialPackingData&packinfo, int p3)
{
//	bool checkoverlap =true;
//	for (int i =0; i<packinfo.numberprocessors; i++)
//	{
//		pthread_t()
//	}
//	return (checkoverlap);
	packinfo.checkoverlap=true;
	double overlapdistance = -packinfo.approximation;
	double distance;
	for (int i=0;i<packinfo.specimen.Particles.Size();i++)
	{
		if (packinfo.particleuses[i]and(!(i==p3)))
		{
			distance = norm(packinfo.specimen.Particles[p3]->x-packinfo.specimen.Particles[i]->x)-packinfo.specimen.Particles[p3]->Dmax-packinfo.specimen.Particles[i]->Dmax;
			if (distance <overlapdistance)
			{
				packinfo.overlappingpoint = i;
				overlapdistance =distance;
				packinfo.checkoverlap=false;
			}
		}
	}
}

inline void Soil::CloseSequentialFace(DEM::SequentialPackingData&packinfo, int face)
{
	packinfo.faces[face].faceuse =0;
	packinfo.numberopenfaces -=1;
}

inline void Soil::CompleteGsdData(DEM::SequentialPackingData&packinfo)
{
    // complete data
    if (packinfo.gsd.intervals[0].begin >0)
    {
    	packinfo.gsd.intervals.insert(packinfo.gsd.intervals.begin(),packinfo.gsd.intervals[0]);
    	packinfo.gsd.intervals[0].begin =0;
    }
    if (packinfo.gsd.intervals[packinfo.gsd.numberintervals-1].begin==100)
    {
    	packinfo.gsd.intervals[packinfo.gsd.numberintervals-2].end=100;
    	packinfo.gsd.intervals[packinfo.gsd.numberintervals-2].enddiameter=packinfo.gsd.intervals[packinfo.gsd.numberintervals-1].begindiameter;
    	packinfo.gsd.intervals.pop_back();
    	packinfo.gsd.numberintervals--;
    }
    else
    {
    	packinfo.gsd.intervals[packinfo.gsd.numberintervals-1].end=100;
    	packinfo.gsd.intervals[packinfo.gsd.numberintervals-1].enddiameter=packinfo.gsd.intervals[packinfo.gsd.numberintervals-1].begindiameter;
    }
    for (int i=0; i<packinfo.gsd.numberintervals-1;i++)
    {
    	packinfo.gsd.intervals[i].end = packinfo.gsd.intervals[i+1].begin;
    	packinfo.gsd.intervals[i].enddiameter=packinfo.gsd.intervals[i+1].begindiameter;
    }
    packinfo.gsd.numberinters =0;
    double count1=0;
    DEM::Interval in;
	for (int i=0; i<packinfo.gsd.numberintervals; i++)
	{
		if ((packinfo.gsd.intervals[i].end -packinfo.gsd.intervals[i].begin)>packinfo.gsd.maxfraction)
		{
			count1 = packinfo.gsd.intervals[i].begin;
			while ((count1 +packinfo.gsd.maxfraction)<packinfo.gsd.intervals[i].end)
			{
				in.begin =count1;
				if (packinfo.gsd.numberinters>0)
				{
					in.begindiameter=packinfo.gsd.inters[packinfo.gsd.numberinters-1].enddiameter;
				}
				else
				{
					in.begindiameter = packinfo.gsd.intervals[0].begindiameter;
				}
				count1 +=packinfo.gsd.maxfraction;
				in.end=count1;
				in.enddiameter= pow(10,log10(packinfo.gsd.intervals[i].begindiameter)+(in.end-packinfo.gsd.intervals[i].begin)/(packinfo.gsd.intervals[i].end-packinfo.gsd.intervals[i].begin)*(log10(packinfo.gsd.intervals[i].enddiameter)-log10(packinfo.gsd.intervals[i].begindiameter)));
				packinfo.gsd.inters.push_back(in);
				packinfo.gsd.numberinters++;
			}
			if (count1 <packinfo.gsd.intervals[i].end)
			{
				in.begin=count1;
				in.begindiameter=packinfo.gsd.inters[packinfo.gsd.numberinters-1].enddiameter;
				count1 = packinfo.gsd.intervals[i].end;
				in.end=count1;
				in.enddiameter=packinfo.gsd.intervals[i].enddiameter;
				packinfo.gsd.inters.push_back(in);
				packinfo.gsd.numberinters++;
			}

		}
		else
		{
			count1=packinfo.gsd.intervals[i].end;
			packinfo.gsd.inters.push_back(packinfo.gsd.intervals[i]);
			packinfo.gsd.numberinters ++;
		}
	}
	if (packinfo.boundarytype=="File")
	{
		// must calculate volume here
	}
	else if (packinfo.boundarytype=="Sphere")
	{
		packinfo.gsd.volume=acos(0.)/3*pow(packinfo.boundarysizes[0],3.);
	}
	else if (packinfo.boundarytype=="Cube")
	{
		packinfo.gsd.volume=pow(packinfo.boundarysizes[0],3.);
	}
	else if (packinfo.boundarytype=="Cylinder")
	{
		packinfo.gsd.volume=acos(0)/4.*pow(packinfo.boundarysizes[0],2.)*packinfo.boundarysizes[1];
	}
	else if (packinfo.boundarytype=="Box")
	{
		packinfo.gsd.volume=packinfo.boundarysizes[0]*packinfo.boundarysizes[1]*packinfo.boundarysizes[2];
	}
	if (packinfo.gsd.specialgravity==0)
	{
		double vol=0.;
		for (int i=0;i<packinfo.gsd.numberinters;i++)
		{
			vol+=(packinfo.gsd.inters[i].end-packinfo.gsd.inters[i].begin)/packinfo.gsd.inters[i].specialgravity;
		}
		packinfo.gsd.mass=100*packinfo.gsd.volume/(vol/(1-packinfo.gsd.porosity));
		packinfo.gsd.density = packinfo.gsd.mass/packinfo.gsd.volume;
		for (int i=0;i<packinfo.gsd.numberinters;i++)
		{
			packinfo.gsd.inters[i].mass= (packinfo.gsd.inters[i].end-packinfo.gsd.inters[i].begin)/100*packinfo.gsd.mass;
			packinfo.gsd.inters[i].volumeparticles=packinfo.gsd.inters[i].mass/packinfo.gsd.inters[i].specialgravity;
		}
	}
	else
	{
		packinfo.gsd.mass=packinfo.gsd.volume*packinfo.gsd.specialgravity*packinfo.gsd.porosity;
		packinfo.gsd.density=packinfo.gsd.mass/packinfo.gsd.volume;
		for (int i=0;i<packinfo.gsd.numberinters;i++)
		{
			packinfo.gsd.inters[i].mass=packinfo.gsd.mass/100*(packinfo.gsd.inters[i].end-packinfo.gsd.inters[i].begin);
			packinfo.gsd.inters[i].specialgravity= packinfo.gsd.specialgravity;
			packinfo.gsd.inters[i].volumeparticles=packinfo.gsd.inters[i].mass/packinfo.gsd.specialgravity;
		}
	}
	for (int i=0;i<10;i++)
	{
		packinfo.boundaryfactors[i]*=packinfo.boundarysizes[i];
	}
	for (int i=0; i<3;i++)
	{
		packinfo.gsd.rectangularboxratios[i]/=packinfo.gsd.rectangularboxratios[0];
	}
}

inline void Soil::CreateSequentialFace( DEM::SequentialPackingData&packinfo, int p[3],int p3, bool sort)
{
	if (sort)
	{
		int count;
		for (int i =0; i<2;i++)
			{
				for (int j=1; j<3; j++)
				{
					if (packinfo.specimen.Particles[p[i]]->Dmax< packinfo.specimen.Particles[p[j]]->Dmax)
						{
							count= p[i];
							p[i]=p[j];
							p[j]=count;
						}
				}
			}
	}
	DEM::SequentialFace facetemprorary;

	Array<Vec3_t> V(3);
	for (int i=0;i<3;i++)
	{
		facetemprorary.points[i]=p[i];
		V[i]= packinfo.specimen.Particles[p[i]]->x;
	}
	// need to calculate constriction ?
	DEM::Face facedemtemprorary(V);															// add face in library
	facedemtemprorary.Normal(facetemprorary.normal);
	facetemprorary.normal /= norm(facetemprorary.normal);
	facetemprorary.point = p3;
	double distance = dot(packinfo.specimen.Particles[p3]->x - packinfo.specimen.Particles[facetemprorary.points[0]]->x, facetemprorary.normal);
	if (distance>0)
	{
		facetemprorary.faceuse = 1;
	}
	else
	{
		facetemprorary.faceuse =-1;
	}
	if (packinfo.particlereuse)
	{
		packinfo.facereuse =false;
		for (int i=0; i< packinfo.numberfaces; i++)
		{
			if (facetemprorary.points[0]== packinfo.faces[i].points[0])
			{
				if ((facetemprorary.points[1]== packinfo.faces[i].points[1])and(facetemprorary.points[2]== packinfo.faces[i].points[2]))
				{
					if (facetemprorary.faceuse*packinfo.faces[i].faceuse <0)
					{
						CloseSequentialFace(packinfo,i);
					}
					else
					{
						packinfo.check =true;
					}
					packinfo.facereuse=true;
					packinfo.facereusenumber = i;
				}
			}
		}
	}
	if (packinfo.facereuse)
	{
		packinfo.facereuse = false;
	}
	else
	{
		packinfo.faces.push_back(facetemprorary);												// add user-defined face
		packinfo.numberfaces +=1;
		packinfo.numberopenfaces +=1;
	}
}

inline void Soil::CreateSequentialTetrahedron(DEM::SequentialPackingData&packinfo, int face, int p3)
{
	DEM::SequentialTetrahedron temprorarytetrahedron;
	temprorarytetrahedron.faces[0]= face;
	int p[3]={packinfo.faces[face].points[0], packinfo.faces[face].points[1], p3};
	CreateSequentialFace(packinfo, p, packinfo.faces[face].points[2]);
	if (packinfo.facereuse)
		{
			temprorarytetrahedron.faces[1]= packinfo.facereusenumber;
			packinfo.facereuse=false;
		}
	else
		{
			temprorarytetrahedron.faces[1]= packinfo.numberfaces-1;
		}
	p[1]= packinfo.faces[face].points[2];
	CreateSequentialFace(packinfo,p, packinfo.faces[face].points[1]);
	if (packinfo.facereuse)
		{
			temprorarytetrahedron.faces[2]= packinfo.facereusenumber;
			packinfo.facereuse=false;
		}
	else
		{
			temprorarytetrahedron.faces[2]= packinfo.numberfaces-1;
		}
	p[0]=packinfo.faces[face].points[1];
	CreateSequentialFace( packinfo,p, packinfo.faces[face].points[0]);
	if (packinfo.facereuse)
		{
			temprorarytetrahedron.faces[3]= packinfo.facereusenumber;
			packinfo.facereuse=false;
		}
	else
		{
			temprorarytetrahedron.faces[3]= packinfo.numberfaces-1;
		}
	packinfo.tetrahedrons.push_back(temprorarytetrahedron);
}

inline void Soil::DrawBoundary( DEM::SequentialPackingData&packinfo)
{
	// Build container
	if (packinfo.boundarytype=="File")
	{
	}
	else if (packinfo.boundarytype=="Sphere")
	{
		packinfo.specimen.AddSphere(-1000,Vec3_t(0.,0.,0.),packinfo.boundaryfactors[0],packinfo.gsd.specialgravity);
	}
	else if (packinfo.boundarytype=="Cube")
	{
		packinfo.specimen.AddCube(-1001,Vec3_t(0.,0.,0.),packinfo.gsd.inters[0].begindiameter/2,packinfo.boundaryfactors[0],packinfo.gsd.specialgravity);
	}
	else if (packinfo.boundarytype=="Cylinder")
	{
		packinfo.specimen.AddCylinder(-1002,Vec3_t(0.,0.,-packinfo.boundaryfactors[1]/2),packinfo.boundaryfactors[0],Vec3_t(0.,0.,packinfo.boundaryfactors[1]/2),packinfo.boundaryfactors[0],packinfo.gsd.inters[0].begindiameter/2,packinfo.gsd.specialgravity);
		packinfo.specimen.AddPlane(-1003,Vec3_t(0.,0.,packinfo.boundaryfactors[1]/2),packinfo.gsd.inters[0].begindiameter/2,1.2*packinfo.boundaryfactors[0],1.2*packinfo.boundaryfactors[0],packinfo.gsd.density,1.);
		packinfo.specimen.AddPlane(-1004,Vec3_t(0.,0.,-packinfo.boundaryfactors[1]/2),packinfo.gsd.inters[0].begindiameter/2,1.2*packinfo.boundaryfactors[0],1.2*packinfo.boundaryfactors[0],packinfo.gsd.density,1.);
	}
	else if (packinfo.boundarytype=="Box")
	{
		packinfo.specimen.GenBox(-1005, packinfo.boundaryfactors[0],packinfo.boundaryfactors[1],packinfo.boundaryfactors[2],packinfo.gsd.inters[0].begindiameter/2,1.2,false);
	}
	for (int i =0; i<packinfo.specimen.Particles.Size(); i++)
	{
		packinfo.boundary.push_back(i);
	}
}

inline void Soil::EstablishLocalSystem( DEM::SequentialPackingData&packinfo, int face)
{
	packinfo.localsystem[0]= packinfo.specimen.Particles[packinfo.faces[face].points[1]]->x - packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x;
	packinfo.localsystem[0]= packinfo.localsystem[0]/norm(packinfo.localsystem[0]);							// e1 vector
	packinfo.localsystem[2]= packinfo.faces[face].normal;													// e3 vector
	packinfo.localsystem[1]= -cross(packinfo.localsystem[0],packinfo.localsystem[2]);
	packinfo.localsystem[1]= packinfo.localsystem[1]/norm(packinfo.localsystem[1]);							// may be not in need
	packinfo.localroot = packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x;
}

inline void Soil::FrozenTag(DEM::SequentialPackingData&packinfo, int tag=-1000)
{
	for (int i=0;i<packinfo.specimen.Particles.Size();i++)
	{
		if (packinfo.specimen.Particles[i]->Tag <= tag)
		{
			packinfo.specimen.Particles[i]->FixVeloc();
		}
	}
}

inline void Soil::PrepareList(DEM::SequentialPackingData&packinfo, bool randomness)
{
	double shapemass;
	double shape;
	double mass;
	double passingmass=0;
	if (randomness)
	{
		for (int i=packinfo.gsd.numberinters-1; i>-1;i--)
		{
			packinfo.gsd.inters[i].numberparticles=0;
			packinfo.gsd.inters[i].mass += passingmass;

		}
	}
	else
	{
		for (int i =packinfo.gsd.numberinters-1; i>-1;i--)
		{
			packinfo.gsd.inters[i].numberparticles=0;
			packinfo.gsd.inters[i].diameter=pow(10,(log10(packinfo.gsd.inters[i].begindiameter)+log10(packinfo.gsd.inters[i].enddiameter))/2);
			packinfo.gsd.inters[i].mass += passingmass;
			shape=0;
			for (int j=0; j<packinfo.gsd.numbershapes;j++)
			{
				shape+=packinfo.gsd.inters[i].shaperatios[j];
			}
			if (shape ==0) throw new Fatal("Error: No proportion for particles in interval from <%s> % to <%d>%", packinfo.gsd.inters[i].begin, packinfo.gsd.inters[i].end);
			mass=0.;
			for (int j=0;j<packinfo.gsd.numbershapes;j++)
			{
				shapemass= ShapeMass(packinfo.gsd, j, packinfo.gsd.inters[i].diameter,packinfo.gsd.inters[i].specialgravity);
				packinfo.gsd.inters[i].shapestates[j][0] = floor(packinfo.gsd.inters[i].shaperatios[j]/shape* (packinfo.gsd.inters[i].mass)/shapemass);
				packinfo.gsd.inters[i].numberparticles += packinfo.gsd.inters[i].shapestates[j][0];
				mass += packinfo.gsd.inters[i].shapestates[j][0]*shapemass;
			}
			passingmass = packinfo.gsd.inters[i].mass -mass;
			packinfo.gsd.inters[i].mass = mass;
			if (packinfo.gsd.inters[i].numberparticles>0)
			{
				packinfo.gsd.inters[i].ability=true;
			}
			else
			{
				packinfo.gsd.inters[i].ability=false;
			}
		}
		// Add particles
		size_t count =0;
		for (int i=packinfo.gsd.numberinters-1; i>-1;i--)
		{
			if (packinfo.gsd.inters[i].ability)
			{
				packinfo.gsd.inters[i].firstparticle=count;
				for (int j=0; j<packinfo.gsd.numbershapes; j++)
				{
					if (packinfo.gsd.inters[i].shapestates[j][0]>0)
					{
						packinfo.gsd.inters[i].shapestates[j][1]=count;
						packinfo.gsd.inters[i].shapestates[j][3]=count;
						for (int k=0; k<packinfo.gsd.inters[i].shapestates[0];k++)
						{
							AddParticle(packinfo, i, j, packinfo.gsd.inters[i].diameter,packinfo.gsd.inters[i].specialgravity,OrthoSys::O,packinfo.gsd.inters[0].begindiameter/2);
							count+=1;
						}
						packinfo.gsd.inters[i].shapestates[j][2]=count-1;
					}
					else
					{
						for (int k=1; k<4;k++)
						{
							packinfo.gsd.inters[i].shapestates[j][k]=-1;
						}
					}
				}
				packinfo.gsd.inters[i].lastparticle=count-1;
				packinfo.gsd.inters[i].usingparticle=packinfo.gsd.inters[i].lastparticle;
			}
		}
		packinfo.gsd.numberparticles = packinfo.specimen.Particles.Size();
	}
}

inline void Soil::PutSequentialParticle(DEM::SequentialPackingData&packinfo, int face, int p3)
{
	EstablishLocalSystem(packinfo,face);															// calculate local coordinates system
	packinfo.checkradius =false;
	double r[4];
	for (int i =0;i<3;i++)
		{
			r[i] = packinfo.specimen.Particles[packinfo.faces[face].points[i]]->Props.R;
		}
	r[3]=packinfo.specimen.Particles[p3]->Props.R;
	double x2 = norm(packinfo.specimen.Particles[packinfo.faces[face].points[1]]->x -packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x);
	double a1 = (r[0] -r[1])/x2;
	double b1 = (pow(x2,2.)+pow(r[0],2.)-pow(r[1],2.))/2./x2;
	double x3 = dot(packinfo.specimen.Particles[packinfo.faces[face].points[2]]->x - packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x,packinfo.localsystem[0]);
	double y3 = dot(packinfo.specimen.Particles[packinfo.faces[face].points[2]]->x - packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x, packinfo.localsystem[1]);
	double a2 = (r[0]-r[2]-a1*x3)/y3;
	double b2 = (pow(x3,2.)+pow(y3,2.)+pow(r[0],2.)-pow(r[2],2.)-2.*b1*x3)/2./y3;
	double x4 = a1*r[3]+b1;
	double y4 = a2*r[3]+b2;
	double z4 = pow(r[0]+r[3],2)-pow(a1*r[3]+b1,2)-pow(a2*r[3]+b2,2);
	if (z4 > 0)
		{
			z4=-packinfo.faces[face].faceuse*pow(z4,0.5);
			Vec3_t position(x4,y4,z4);
			Vec3_t finalposition =  packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x;
			for (int i=0; i<3; i++)
				{
					for (int j=0; j<3; j++)
						{
							finalposition(i)+= position(j)*packinfo.localsystem[j](i);
						}
				}
			packinfo.specimen.Particles[p3]->Position(finalposition);
			packinfo.checkradius =true;
		}
}

inline void Soil::ReadGsd( DEM::SequentialPackingData&packinfo, string FileKey)
{
	string filename=FileKey;
	filename.append(".gsd");
	if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.c_str());
    ifstream datain;
    datain.open(filename.c_str());
    printf("\n%s--- Loading GSD file %s --------------------------------------------%s\n",TERM_CLR1,filename.c_str(),TERM_RST);
    datain >> packinfo.gsd.Soilname;					datain.ignore(200,'\n');
    datain >> packinfo.boundarytype;						datain.ignore(200,'\n');
    if (packinfo.boundarytype =="File")
    {
    	datain >>packinfo.gsd.boundaryfile;			datain.ignore(200,'\n');
    	packinfo.gsd.boundaryfile.append(".h5");
    	if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",packinfo.gsd.boundaryfile.c_str());
    	// how to check overlap with boundary object ?
    	packinfo.specimen.Load(packinfo.gsd.boundaryfile.c_str());
    }
    else
    {
        int numbersize;
        datain >> numbersize;
    	for (int i=0;i<numbersize;i++)
    	{
    	   	datain >> packinfo.boundarysizes[i];
    	}
    	for (int i=0; i<numbersize;i++)
    	{
    		datain >> packinfo.boundaryfactors[i];
    	}
    }
    										datain.ignore(200,'\n');
    datain >> packinfo.gsd.arrangement;				datain.ignore(200,'\n');
    datain >> packinfo.numberprocessors; 	datain.ignore(200,'\n');
    datain >> packinfo.gsd.specialgravity;			datain.ignore(200,'\n');
    datain >> packinfo.gsd.porosity;					datain.ignore(200,'\n');
    datain >> packinfo.gsd.maxfraction;				datain.ignore(200,'\n');
    datain >> packinfo.gsd.numberintervals;
    datain >> packinfo.gsd.numbershapes;				datain.ignore(200,'\n');
    cout << packinfo.gsd.numbershapes <<"\n";
    DEM::Interval inter;
    for (int i=0; i<packinfo.gsd.numberintervals; i++)
    {
    	datain >> inter.begin;
    	datain >> inter.begindiameter;
    	datain >> inter.specialgravity;
    	for (int j=0; j <packinfo.gsd.numbershapes; j++)
    	{
    		datain >> inter.shaperatios[j];
    		cout << inter.shaperatios[j];
    	}
    										datain.ignore(200,'\n');
    	packinfo.gsd.intervals.push_back(inter);
    }
    datain >> packinfo.approximation;
    for (int i =0; i<3; i++)
    {
    	datain >> packinfo.gsd.rectangularboxratios[i];
    }
    										datain.ignore(200,'\n');
    datain.close();
}

inline void Soil::SequentialPacking(DEM::SequentialPackingData&packinfo)
{
//	for (int i=0; i<packinfo.numberprocessors; i++)
//	{
//		packinfo.threadnumbers.push_back(int(i/packinfo.numberprocessors*specimen.Particles.Size()));
//	}
//	packinfo.threadnumbers.push_back(specimen.Particles.Size());
	packinfo.particlereuse =false;
	packinfo.facereuse =false;
	packinfo.check=true;
	for (int i=0; i<packinfo.gsd.numberparticles; i++)
	{
		packinfo.particleuses.push_back(false);
	}
	if (packinfo.gsd.arrangement=="Layer")
	{

	}
	else if (packinfo.gsd.arrangement=="Discrete")
	{
		packinfo.usinginter = packinfo.gsd.numberinters-1;
		packinfo.usingparticle = packinfo.specimen.Particles.Size()-1;
		int p[4]= {packinfo.usingparticle, packinfo.usingparticle-1, packinfo.usingparticle-2, packinfo.usingparticle-3};
		BasicTetrahedron(packinfo,p);	// create initial tetrahedron

		packinfo.usingparticle -=4;
		// packing loop
		while ((packinfo.numberopenfaces >0)and(packinfo.numberunusedparticles >0))
		{
			while (packinfo.particleuses[packinfo.usingparticle])
				{
					packinfo.usingparticle -=1;
					if (packinfo.usingparticle==0)
						{
							break;
						}
					}
			while (packinfo.faces[packinfo.usingface].faceuse == 0)
				{
					packinfo.usingface +=1;
					if (packinfo.usingface == packinfo.numberfaces)
						{
							break;
						}
				}
			TrySequentialParticle( packinfo, packinfo.usingface, packinfo.usingparticle);
		}
	}
}

inline void Soil::TrySequentialParticle( DEM::SequentialPackingData&packinfo, int face, int p3)
{
	packinfo.temproraryparticle = p3;
	double distancepoint =0;
	double distanceoverlap =0;
	PutSequentialParticle(packinfo,face,p3);
	CheckOverlap( packinfo, p3);
	CheckBoundary(packinfo,p3);
	if ((!packinfo.checkradius)or (!packinfo.checkboundary) or (!packinfo.checkoverlap))
	{
		int interval = floor(abs(packinfo.specimen.Particles[packinfo.temproraryparticle]->Tag)/packinfo.gsd.numbershapes);
	}
}

inline void Soil::UseParticle( DEM::SequentialPackingData&packinfo, int p3)
{
	if (packinfo.particleuses[p3])
	{
		cout << "Particle " << p3 << " has been used before \n";
	}
	packinfo.particleuses[p3]=true;
	packinfo.numberunusedparticles-=1;
	int inter = floor(abs(packinfo.specimen.Particles[p3]->Tag/packinfo.gsd.numbershapes));
	packinfo.gsd.inters[inter].usingparticle -=1;
	if (packinfo.gsd.inters[inter].usingparticle < packinfo.gsd.inters[inter].firstparticle )
		{
			packinfo.gsd.inters[inter].ability =false;
			if (inter == packinfo.usinginter)
				{
					packinfo.usinginter -=1;
				}
		}
}
}
int main(int argc, char **argv)
{
	DEM::Soil sand;
	DEM::SequentialPackingData packinfo;
	sand.ReadGsd(packinfo);
	sand.CompleteGsdData(packinfo);
	sand.PrepareList(packinfo, false);
	sand.SequentialPacking(packinfo);
}
