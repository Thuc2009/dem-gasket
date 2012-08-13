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
    // domain
    DEM::Domain dom;

    String filename(argv[1]);

    dom.Load(filename.CStr());


    //Erase the plane
    Array<int> DeleteTags(1);
    DeleteTags = -1;   // These tags are of all the obstacles
    dom.DelParticles(DeleteTags);

    //Add the container
    dom.GenBoundingBox (-1,0.1,1.2,false);

    //Fixing the plates so they dont move
    for (size_t i=-1;i>-7;i--)
    {
        dom.GetParticle(i)->FixVeloc();
    }



    // Set the parameters
    double Kn          = 1.0e6;
    double Kt          = 3.3e5;
    double Gn          = 16.0;
    double Gt          = 8.0;
    double Mu          = 0.2;

    Dict P;
    for (size_t i=0;i<20;i++)
    {
        P.Set(/*Tag*/i,"Kn Kt Gn Gt Mu",Kn,Kt,Gn,Gt,Mu);    
    }
    dom.SetProps(P);

    //Add the gravity
    Vec3_t g(0.0,0.0,-9.8);
    for (size_t i=0;i<dom.Particles.Size();i++)
    {
        dom.Particles[i]->Ff = dom.Particles[i]->Props.m*g;
    }

    dom.Solve     (/*tf*/10.0, /*dt*/1.0e-3, /*dtOut*/0.1, NULL, NULL, /*filekey*/"dem",/*Visit visualization*/3,/*N_proc*/3);

    return 0;
}
MECHSYS_CATCH
