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
#include <mechsys/dem/domain.h>
// variable declaration
using namespace std;
using namespace DEM;
// declare using structure
// global data
struct interval																// structure of interval of GSD
{
	bool ability;
	int begin;																// fraction
	int end;
	int firstparticle;														// particle number
	int lastparticle;
	int usingparticle;
	double density;
	double begindiameter;
	double enddiameter;
	double diameter;														//mean diameter;
	double radius;
	double mass;
	int numberparticles;
	int shapeparticles[10];
	int unusedshapeparticles[10];
	double shapefactors[10];
	double shapemass[10];
	double equalradius[10];
	double equaldiameter[10];
	double rectangularboxsize[3];
};
struct data
{
	bool choice;
	bool check;
	bool facereuse;
	bool particlereuse;
	bool checkoverlap;
	bool checkdistance;
	bool checkradius;
	char choicechar;
	double approximation;
	double density;																			// soil density might not be a constant
	double mass;
	double porosity;
	double volume;
	int facereusenumber;
	int count;
	int maxfraction;
	int GSD;
	int numberintervals;
	int numberparticles;
	int numberunusedparticles;
	int numberfaces;
	int numberopenfaces;
	int numbershapes;
	int overlappingpoint;
	int usingface;
	int usinginterval;
	int usingparticle;
	string datafile;
	string outputdomain;
	string soilname;
	string line;
	string specimentype;
	string fillingtype;
	string particleshapetype;
	double specimen[3];
	double factor[3];
	double rectangularboxratio[3];
	double gsd[101];
	double densities[101];
	double shapefactors[101][10];
	double csd_number[101];
	double csd_area[101];
	double diameters[101];
	Vec3_t localroot;
	Vec3_t localsystem[3];
	Vec3_t position;
	Vec3_t finalposition;
	vector <int> gsd1;
	DEM::Domain particles;
	vector<interval> intervals;
	vector<face> faces;
	//vector<DEM::Face> demfaces;
	vector<bool> particlesuse;
	vector<tetrahedron> tetrahedrons;
};
// function and procedure declaration
void calculateintervals(data & dat);
void checkoverlap(data&dat, int p0);
void checkdistance(data&dat, int p0);
void closeface(data&dat, int m0);
void convertcoordinates(data & dat, Vec3_t & position, bool way);
void createface(data & dat, int p0, int p1, int p2,int faceuse);
void createtetrahedron(data & dat, int m0, int p0);
void establishlocalsystem(data&dat, int m0);
void gsdgenerate(data & dat);
void listprepare(data & dat);
void moveon(data & dat, int p0);
void putparticle(data & dat, int p0, int m0);
void savedomain(data&dat);
void textout (data& dat);
void useparticle(data & dat, int p0);
double equalradius(data&dat, int p0);
double interpolation (double x1,double x2, double y1,double y2,double x3);
double getinterval (data&dat, int p0);
// main procedure
int main(int argc, char **argv) try
	{
	data dat;
	dat.check = true;
	dat.particlereuse = false;
	dat.facereuse = false;
	dat.approximation =0.001;
	int fraction=0;
	int interval=0;
	int usingparticle =0;
	int usingface = 0;
	double distance1=0;
	double distance2=0;
	ifstream datain;
	//cout << "Welcome to CSD program \n";
	//cout << "Enter data file name: ";
	//cin >> dat.datafile;
	dat.datafile = "test.csd";
	datain.open(dat.datafile.c_str());
	//if (!datain.good())
	//	{
	//		dat.datafile = "sand.csd";
	//		datain.open(dat.datafile.c_str());
	//	}
	datain >> dat.soilname; 		datain.ignore(200,'\n');
	datain >> dat.specimentype;		datain.ignore(200,'\n');
	datain >> dat.particleshapetype;datain.ignore(200,'\n');
	datain >> dat.fillingtype;		datain.ignore(200,'\n');
	datain >> dat.density; 			datain.ignore(200,'\n'); 																								// density of soil might not a constant
	datain >> dat.porosity;			datain.ignore(200,'\n');
	datain >> dat.maxfraction; 		datain.ignore(200,'\n');
	datain >> dat.volume;			datain.ignore(200,'\n');
	datain >> dat.GSD;
	datain >> dat.numbershapes;		datain.ignore(200,'\n');
	for (int i=1; i<dat.GSD+1;i++)
		{
			datain >> fraction;
			datain >> dat.gsd[fraction];
			dat.gsd1.push_back(fraction);
			datain >> dat.densities[fraction];
			for (int j =0; j<dat.numbershapes;j++)
				{
					datain >>dat.shapefactors[fraction][j];
				}
			datain.ignore (200,'\n');
		}
	for (int i=0; i<3; i++)
		{
			datain >> dat.factor[i];
		}
									datain.ignore (200,'\n');
	datain >> dat.outputdomain;		datain.ignore (200,'\n');
	for (int i=0; i<3; i++)
		{
			datain >> dat.rectangularboxratio[i];
		}
									datain.ignore (200,'\n');
	datain.close();
	gsdgenerate(dat);																																		// Generate gsd
	// choose again volume
	dat.choice = false;
	while (!dat.choice)
		{
			calculateintervals(dat);																														// generate list of intervals with number of particles
			// check information after calculation
			for (int i =0; i<dat.numberintervals;i++)
				{
					cout << dat.intervals[i].begin << " " << dat.intervals[i].end<< " "<< dat.intervals[i].density<< " "<< dat.intervals[i].diameter<< " ***** ";
					for (int j=0; j<dat.numbershapes; j++)
						{
							cout << dat.intervals[i].shapefactors[j]<<" ";
						}
					cout << " ******* ";
					for (int j=0; j<dat.numbershapes; j++)
						{
							cout << dat.intervals[i].shapemass[j]<<" ";
						}
					cout << " ******* ";
					for (int j=0; j<dat.numbershapes; j++)
						{
							cout << dat.intervals[i].shapeparticles[j]<<" ";
						}
					cout << dat.intervals[i].mass << "\n";
				}
			cout << "mass: " << dat.mass<<"\n";
			for (int i=0; i<dat.numberintervals; i++)
				{
					cout << i << " " << dat.intervals[i].begin<< " " << dat.intervals[i].numberparticles << " " << dat.intervals[i].diameter << "\n";
				}
			cout << dat.numberparticles;
			//cout << "Do you want to use this input ? Y/n: ";
			//cin >> dat.choicechar;
			//if (dat.choicechar == 'Y')
			//	{
			//		cout << dat.choicechar;
			//		dat.choice=true;
			//	}
			//else
			//	{
			//		cout << "Please enter new value of volume: ";
			//		cin >> dat.volume;
			//	}
			dat.choice =true; 																	// autocheck
		}
	//generate specimen;
	if (dat.specimentype.compare("Cube")==0)
		{
			for (int i =0; i<3; i++)
				{
					dat.specimen[i]= pow(dat.volume,1.0/3.0)/2;
					dat.specimen[i]*=dat.factor[i];
				}

			//Vec3_t dummy(dat.specimen[0],0,0);
			//dat.particles.AddPlane(-1,dummy,0.1,2*dat.specimen[0],2*dat.specimen[0],1.0,0.5*M_PI,&OrthoSys::e1);
		}
	else if (dat.specimentype.compare("Sphere")==0)
		{
			for (int i = 0; i<3; i++)
				{
					dat.specimen[i]=pow((3.0 * dat.volume /8.0/ acos(0)),(1.0/3.0));
					dat.specimen[i]*=dat.factor[i];
				}
		}
	dat.numberunusedparticles = dat.numberparticles;
	// Generate particles

	for (int i=0; i< dat.particles.Particles.Size(); i++)
		{
			if (dat.particlesuse[i]==false)
				{
					dat.particles.Particles[i]->Tag -=100*dat.numbershapes;
				}
		}
	// export results
	savedomain(dat);
	textout(dat);
	cout << "saved domain \n";
    return 0;
}
MECHSYS_CATCH
// function
inline double equalradius(data&dat, int p0)
	{
		double r= dat.intervals[abs(int(dat.particles.Particles[p0]->Tag/dat.numbershapes))].equalradius[abs(dat.particles.Particles[p0]->Tag)%dat.numbershapes];
		return (r);
	}
inline double interpolation (double x1,double x2, double y1,double y2,double x3)
	{
		double y3;
		y3= y1 + (x3-x1)*(y2-y1)/(x2-x1);
		return (y3);
	}
inline double getinterval(data&dat, int p0)
	{
		double interval = abs(int(dat.particles.Particles[p0]->Tag/dat.numbershapes));
		return (interval);
	}
// procedure
inline void calculateintervals(data & dat)																							//prepare a list of intervals
	{
		// generate intervals
		interval inter;
		dat.numberintervals =0;
		int count = 0;
		if (dat.gsd1[0]>0)
			{
				while (count < dat.gsd1[0])
					{
						dat.numberintervals ++;
						inter.begin =count;
						dat.intervals.push_back(inter);
						count += dat.maxfraction;
					}
			}
		for (int i=0; i<dat.GSD; i++)
			{
				if (dat.maxfraction>(dat.gsd1[i+1]-dat.gsd1[i]))
					{
						dat.numberintervals ++;
						inter.begin = dat.gsd1[i];
						dat.intervals.push_back(inter);
					}
				else
					{
						count = dat.gsd1[i];
						while (count < dat.gsd1[i+1])
							{
								dat.numberintervals ++;
								inter.begin = count;
								dat.intervals.push_back(inter);
								count += dat.maxfraction;
							}
					}
			}
		if (dat.intervals[dat.numberintervals-1].begin =100)
			{
				dat.intervals.pop_back();																												// not completely delete, therefore, still active below
				dat.numberintervals --;
			}
		dat.mass = dat.density * dat.volume * (1-dat.porosity);																							// if density of soil is constant
		double module =0.0;																																// abundant mass
		dat.numberparticles =0;
		for (int i =dat.numberintervals-1; i> -1; i--)
			{
				dat.intervals[i].end = dat.intervals[i+1].begin;																						// not completely delete before
				dat.intervals[i].begindiameter = dat.gsd[dat.intervals[i].begin];
				dat.intervals[i].enddiameter = dat.gsd[dat.intervals[i].end];
				dat.intervals[i].density= dat.densities[dat.intervals[i].begin]; 																									// if density of soil is constant
				dat.intervals[i].mass = dat.mass*(dat.intervals[i].end-dat.intervals[i].begin)/100.0 + module;
				dat.intervals[i].diameter = pow(10.0, (log10(dat.intervals[i].begindiameter)+log10(dat.intervals[i].enddiameter))/2); 					// generate mean diameter of interval
				dat.intervals[i].radius = dat.intervals[i].diameter/2.0;
				dat.intervals[i].rectangularboxsize[2] = pow(2.0, 0.5)*dat.intervals[i].diameter*dat.rectangularboxratio[2]/(dat.rectangularboxratio[1]+dat.rectangularboxratio[2]);
				dat.intervals[i].rectangularboxsize[0] =dat.rectangularboxratio[0]/dat.rectangularboxratio[2]*dat.intervals[i].rectangularboxsize[2];
				dat.intervals[i].rectangularboxsize[1] =dat.rectangularboxratio[1]/dat.rectangularboxratio[2]*dat.intervals[i].rectangularboxsize[2];
				dat.intervals[i].shapemass[0] = dat.intervals[i].density * pow(dat.intervals[i].diameter,3.0)/3.0*acos(0.0);												// mass of one particle must be transfer to g/cm3
				dat.intervals[i].shapemass[1] = dat.intervals[i].density*pow(dat.intervals[i].diameter,3.0);																// cube
				dat.intervals[i].shapemass[2] = dat.intervals[i].density*pow(dat.intervals[i].diameter,3.0)/3.0;															// tetrahedron, the edge is the diagonal of net
				dat.intervals[i].shapemass[3] = dat.intervals[i].density*dat.intervals[i].rectangularboxsize[0]*dat.intervals[i].rectangularboxsize[1]*dat.intervals[i].rectangularboxsize[2];  																				// Rectangular box
				double factorsum =0.0;
				for (int j =0; j<dat.numbershapes; j++)
					{
						dat.intervals[i].shapefactors[j]=dat.shapefactors[dat.intervals[i].begin][j];
						factorsum += dat.intervals[i].shapefactors[j];
					}
				for (int j=0; j<dat.numbershapes; j++)
					{
						dat.intervals[i].shapefactors[j]/=factorsum;
						dat.intervals[i].shapeparticles[j]=dat.intervals[i].mass*dat.intervals[i].shapefactors[j]/dat.intervals[i].shapemass[j];							// number of shaped particles
						dat.intervals[i].unusedshapeparticles[j]=dat.intervals[i].shapeparticles[j];													// number of unused shaped particles
					}
				dat.intervals[i].equalradius[0]=dat.intervals[i].radius;
				dat.intervals[i].equalradius[1]=dat.intervals[i].radius*pow(3.0,0.5);
				dat.intervals[i].equalradius[2]=dat.intervals[i].radius*pow(3.0,0.5);
				dat.intervals[i].equalradius[3]=norm(Vec3_t(dat.intervals[i].rectangularboxsize[0],dat.intervals[i].rectangularboxsize[1],dat.intervals[i].rectangularboxsize[2]))/2;
				dat.intervals[i].numberparticles =0;
				module = dat.intervals[i].mass;
				for (int j=0; j<dat.numbershapes;j++)
					{
						dat.intervals[i].equaldiameter[j]=2.0*dat.intervals[i].equalradius[j];
						dat.intervals[i].numberparticles += dat.intervals[i].shapeparticles[j];															// number of particles in interval
						module -= dat.intervals[i].shapeparticles[j]*dat.intervals[i].shapemass[j];																		// abundant mass is passed to next interval
					}
				dat.intervals[i].mass = dat.intervals[i].mass -module;																					// real mass of interval after calculation
				dat.numberparticles += dat.intervals[i].numberparticles;
				if (dat.intervals[i].numberparticles > 0)
					{
						dat.intervals[i].ability=true;
					}
				else
					{
						dat.intervals[i].ability =false;
					}
			}
	}
inline void gsdgenerate(data & dat)
	{
		for (int i =0; i<dat.GSD; i++)
				{
					for (int j=dat.gsd1[i]+1; j<dat.gsd1[i+1];j++)
						{
							dat.gsd[j] = interpolation (dat.gsd1[i],dat.gsd1[i+1],dat.gsd[dat.gsd1[i]],dat.gsd[dat.gsd1[i+1]],j);  						//diameter might be generated by log function
							dat.densities[j] = dat.densities[dat.gsd1[i]];
							for (int k =0; k<dat.numbershapes; k++)
								{
									dat.shapefactors[j][k]=dat.shapefactors[dat.gsd1[i]][k];
								}
						}
				}
	}
inline void listprepare(data & dat)																								// add particles to domain
	{
		int count =0;
		int count1 =0;
		Vec3_t x(0.,0.,0.);
		for (int i=0;i<dat.numberintervals;i++)
			{
				if (dat.intervals[i].ability == true)
					{
						dat.intervals[i].firstparticle = count;
						for (int j=0; j<dat.numbershapes; j++)
							{
								for (int k=0; k < dat.intervals[i].shapeparticles[j]; k++)
									{
										if (j==0)
											{
												dat.particles.AddSphere(i*dat.numbershapes, x ,dat.intervals[i].radius, dat.intervals[i].density); 				//x - coordinates has been declared before
											}
										else if (j==1)
											{
												dat.particles.AddCube(-i*dat.numbershapes-j,x,0.05*dat.intervals[i].radius,dat.intervals[i].radius,dat.intervals[i].density);
											}
										else if (j==2)
											{
												dat.particles.AddTetra(-i*dat.numbershapes-j,x,0.05*dat.intervals[i].radius,pow(2.0,0.5)*dat.intervals[i].diameter, dat.intervals[i].density);
											}
										else if (j==3)
											{
												dat.particles.AddRecBox(-i*dat.numbershapes -j, x,Vec3_t(dat.intervals[i].rectangularboxsize[0],dat.intervals[i].rectangularboxsize[1],dat.intervals[i].rectangularboxsize[2]),0.05*dat.intervals[i].radius,dat.intervals[i].density);
											}
										dat.particlesuse.push_back(false);
										count +=1;
									}
							}
						dat.intervals[i].lastparticle =count -1;
						dat.intervals[i].usingparticle =count-1;
						dat.intervals[i].ability =true;
					}
			}
	}
inline void savedomain(data&dat)
	{
		dat.particles.WriteXDMF(dat.outputdomain.c_str());// export to draw visual results
		dat.particles.Save(dat.outputdomain.c_str());
	}
inline void textout (data & dat)
	{
		ofstream dataout;																							// export to text file
		dataout.open("sand.out");
		dataout << dat.numberparticles << " number of particles" << "\n";
		dataout << dat.numberunusedparticles << " number of unused particles" << "\n";
		dataout << dat.numberfaces << " number of faces" << "\n";
		dataout << dat.numberopenfaces << " number of opening faces" << "\n";
		dataout << dat.specimentype << " " << dat.volume << " " << dat.specimen[0] << "\n";
		for (int i=0; i<dat.numberintervals; i++)
			{
				dataout << i << " " << dat.intervals[i].diameter<< " " << dat.intervals[i].numberparticles <<" "<< dat.intervals[i].mass << " " << dat.intervals[i].firstparticle << " "<< dat.intervals[i].lastparticle<< "\n";
			}
		dataout << "List of unused particles \n";
		for (int i=0; i<dat.numberintervals; i++)
			{
				if (dat.intervals[i].ability == true)
					{
						dataout << i << " : " << dat.intervals[i].usingparticle << " " <<dat.intervals[i].firstparticle << "\n";
					}
			}
		dataout.close();
	}
