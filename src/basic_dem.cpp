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
    dom.AddCube(-1,Vec3_t(5,0.0,5.0),0.1,1.0,3.0,M_PI/3.0,&OrthoSys::e0);
    dom.AddTetra(-2,Vec3_t(-5,0.0,5.0),0.1,1.5,3.0,M_PI/3.0,&OrthoSys::e0);

    dom.AddPlane(-3,OrthoSys::O,0.1,20.0,20.0,3.0);
    dom.GetParticle(-3)->FixVeloc();

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
