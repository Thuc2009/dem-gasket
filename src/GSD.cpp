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
#include </home/thuc2009/dem-gasket/src/gsd.h>
// using namespace
using namespace std;
using namespace DEM;
// this namespace


int main(int argc, char **argv)
{
	srand(time(NULL));
//	int s=rand()%100;
//	cout << s<<endl;
	DEM::Soil sand;
	DEM::SequentialPackingData packinfo;
	string path;
	if (argc<2)
		{
			path = "./1/1.gsd";
		}
	else
		{
			path = argv[1];
		}
	sand.ReadGsd(packinfo, path);
	cout << "Start Time: "<<sand.Now() << endl;
	sand.CompleteGsdData(packinfo);
	sand.FindMinimumSpecimenSize(packinfo);
	sand.PrepareList(packinfo);
	cout << "Number of particles: "<<packinfo.gsd.numberparticles<<endl;
	// sand.PrintOut(packinfo);
	sand.SequentialPacking(packinfo);
	sand.DeleteUnusedParticles(packinfo);
	cout << "End time: "<<sand.Now()<<endl;
	//packinfo.specimen.AddCylinder(-1000,Vec3_t(0.,0.,-packinfo.boundaryfactors[1]/2),packinfo.boundaryfactors[0]/2,Vec3_t(0.,0.,packinfo.boundaryfactors[1]/2),packinfo.boundaryfactors[0]/2,0.1,2.65);
	sand.SaveDomain(packinfo.specimen,packinfo.gsd.Soilname,1);
	sand.TextOut(packinfo);

//	sand.ReadDemParameters(packinfo,"sergio");
//	sand.DropDown(packinfo);
//	sand.SaveDomain(packinfo.specimen,packinfo.gsd.Soilname,1);
}
