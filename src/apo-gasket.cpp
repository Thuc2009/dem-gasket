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

void ApoGasket(DEM::Domain & d, double L)
{

    // Put the first particle
    d.Particles[0]->Position(Vec3_t(0.0,0.0,0.0));

    //Put the second particle
    d.Particles[1]->Position(Vec3_t(0.5*L,0.5*L,0.5*L));

    //Put the third particle
    d.Particles[2]->Position(Vec3_t(-0.5*L,-0.5*L,-0.5*L));
}


int main(int argc, char **argv) try
{
    // domain
    DEM::Domain dom;

    Vec3_t x1(0.0,0.0,0.0);
    dom.AddSphere(-1,x1/*position*/,1.0/*radius*/,3.0/*density*/);
    Vec3_t x2(10.0,0.0,0.0);
    dom.AddSphere(-2,x2/*position*/,1.0/*radius*/,3.0/*density*/);
    Vec3_t x3(10.0,10.0,0.0);
    dom.AddSphere(-3,x3/*position*/,1.0/*radius*/,3.0/*density*/);

    ApoGasket(dom,10.0);

    dom.WriteXDMF("apo-gasket");

    return 0;
}
MECHSYS_CATCH
