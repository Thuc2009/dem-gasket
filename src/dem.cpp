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


int main(int argc, char **argv) try
{
    if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);

    size_t Nproc = 1; 
    if (argc==3) Nproc=atoi(argv[2]);
    String filekey  (argv[1]);
    String filename (filekey+".par");
    String fileout;
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    ifstream infile(filename.CStr());
    // Set the parameters

    double Alpha;
    double Kn;
    double Kt;
    double Gn;
    double Gt;
    double Mu;
    double Press;
    double dt;
    double dtOut;
    double Tf;

    infile >> filename;     infile.ignore(200,'\n');
    infile >> fileout;      infile.ignore(200,'\n');
    infile >> Alpha;        infile.ignore(200,'\n');
    infile >> Kn;           infile.ignore(200,'\n');
    infile >> Kt;           infile.ignore(200,'\n');
    infile >> Gn;           infile.ignore(200,'\n');
    infile >> Gt;           infile.ignore(200,'\n');
    infile >> Mu;           infile.ignore(200,'\n');
    infile >> Press;        infile.ignore(200,'\n');
    infile >> dt;           infile.ignore(200,'\n');
    infile >> dtOut;        infile.ignore(200,'\n');
    infile >> Tf;           infile.ignore(200,'\n');

    // domain
    DEM::Domain dom;

    dom.Load(filename.CStr());

    dom.GenBoundingBox (/*InitialTag*/-1000, 0.05, /*Cf*/1.1);
    dom.GetParticle(-1000)->FixVeloc();
    dom.GetParticle(-1001)->FixVeloc();
    dom.GetParticle(-1002)->FixVeloc();
    dom.GetParticle(-1003)->FixVeloc();
    dom.GetParticle(-1004)->FixVeloc();
    dom.GetParticle(-1005)->FixVeloc();

    double area = (dom.GetParticle(-1000)->x(0) - dom.GetParticle(-1001)->x(0))*(dom.GetParticle(-1002)->x(1) - dom.GetParticle(-1003)->x(1));



    dom.GetParticle(-1004)->vzf = false;
    dom.GetParticle(-1005)->vzf = false;
    dom.GetParticle(-1004)->Ff  = Vec3_t(0.0,0.0,-Press*area);
    dom.GetParticle(-1005)->Ff  = Vec3_t(0.0,0.0, Press*area);


    //dom.WriteXDMF("dem");


//
    //Add the gravity
    //Vec3_t g(0.0,0.0,9.8);
    //

    double minmass = 10.0e12;
    double minsr   = 10.0e12;
    double mindiam = 10.0e12;
    
    for (size_t i=0;i<dom.Particles.Size();i++)
    {
        dom.Particles[i]->Props.Kn = Kn;
        dom.Particles[i]->Props.Kt = Kt;
        dom.Particles[i]->Props.Gn = Gn;
        dom.Particles[i]->Props.Gt = Gt;
        dom.Particles[i]->Props.Mu = Mu;
        if (dom.Particles[i]->Dmax    < mindiam) mindiam = dom.Particles[i]->Dmax   ;
        if (dom.Particles[i]->Props.R < minsr  ) minsr   = dom.Particles[i]->Props.R;
        if (dom.Particles[i]->Props.m < minmass) minmass = dom.Particles[i]->Props.m;
    }

    std::cout << "Minimun size             " << mindiam          << std::endl;
    std::cout << "Minimun sphero-radius    " << minsr            << std::endl;
    std::cout << "Minimun mass             " << minmass          << std::endl;
    std::cout << "Suggested time step      " << sqrt(minmass/Kn) << std::endl;

//
    dom.Solve     (/*tf*/Tf, /*dt*/dt, /*dtOut*/dtOut, NULL, NULL, /*filekey*/filekey.CStr(),/*Visit visualization*/2,/*N_proc*/Nproc);
    dom.Save(fileout.CStr());

    return 0;
}
MECHSYS_CATCH
