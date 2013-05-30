/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2012 To Huu Duc, Sergio Torres                         *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

// MechSys
#include <mechsys/dem/domain.h>
using namespace std;
using namespace DEM;
void Report (DEM::Domain & dom, void *UD)
{
    dom.Save(dom.FileKey.c_str());
    //std::cout << dom.Time << std::endl;
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
	double dtOut;
	double Alpha;
	double Kinematicenergy;
	int starttag;
	int numbershapes;
	int numberintervals;
	string filekey;
	int visualization;
	int processornumber;
	double roundratio;
	DEM::Domain dom;
	//DEM::Domain particlesbefore;
	//Array<Vec3_t> posbefore;
	string datafile="suffusion.drp";
	string domainin;
	string domainout;
	ifstream datain;
	datain.open(datafile.c_str());
	datain >> domainin; 				datain.ignore(200,'\n');
	datain >> domainout;				datain.ignore(200,'\n');
	datain.close();
	dom.Load(domainin.c_str());
	datain.open("suffusion.par");
	datain >> Kn;				datain.ignore(200,'\n');
	datain >> Kt;				datain.ignore(200,'\n');
	datain >> Gn;				datain.ignore(200,'\n');
	datain >> Gt;				datain.ignore(200,'\n');
	datain >> Mu;				datain.ignore(200,'\n');
	datain >> dt;				datain.ignore(200,'\n');
	datain >> tf;				datain.ignore(200,'\n');
	datain >> dtOut;			datain.ignore(200,'\n');
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
    // domain
    Vec3_t g(0.0,9.8,-9.8);
    Vec3_t xmin,xmax;
    dom.BoundingBox(xmin,xmax);
    dom.GenBox(-1000, 1.0*(xmax(0)-xmin(0))+0.1,2.0*(xmax(1)-xmin(1))+0.1, 1.0*(xmax(2)-xmin(2))+0.1, 0.05, 1.1);
    dom.GetParticle(-1000)->FixVeloc();
    dom.GetParticle(-1001)->FixVeloc();
    dom.GetParticle(-1002)->FixVeloc();
    dom.GetParticle(-1003)->Position(Vec3_t(0,xmin(1)-0.05,0));
    dom.GetParticle(-1003)->FixVeloc();
    dom.GetParticle(-1004)->FixVeloc();
    dom.GetParticle(-1005)->FixVeloc();
    for (size_t i=0;i<dom.Particles.Size();i++)
    {
        if (dom.Particles[i]->Tag>-1000)
    	{
            dom.Particles[i]->Props.Kn = Kn;
            dom.Particles[i]->Props.Kt = Kt;
            dom.Particles[i]->Props.Gn = Gn;
            dom.Particles[i]->Props.Gt = Gt;
            dom.Particles[i]->Props.Gv = 0.01;
            dom.Particles[i]->Props.Mu = Mu;
        	if (dom.Particles[i]->Props.R > 0.1)
        	{
        		dom.Particles[i]->FixVeloc();
        	}
        	else
        	{
                dom.Particles[i]->Ff = dom.Particles[i]->Props.m*g;
        	}
    	}
    }
    dom.Solve(/*tf*/tf, /*dt*/dt, /*dtOut*/dtOut, NULL, NULL, /*filekey*/filekey.c_str(),/*Visit visualization*/visualization,/*N_proc*/processornumber);
    dom.Save(domainout.c_str());
    dom.WriteXDMF(domainout.c_str());
    return 0;
}
MECHSYS_CATCH
