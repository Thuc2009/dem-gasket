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
struct tetrahedron
{
	int points[4];
	int faces[4];
	double poreradius;
	Vec3_t porecentre;
};
struct face
{
	int points[3];
	int point;
	int faceuse;
	double constrictionsize;
	Vec3_t normal;
	Vec3_t constriction;
	Vec3_t constrictioncentre;
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
	int temproraryparticle;
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
void moveon(data & dat, int p0, int m0);
void putparticle(data & dat, int p0, int m0);
void savedomain(data&dat);
void textout (data& dat);
void tryputparticle(data&dat, int p0,int m0);
void useparticle(data & dat, int p0);
double equalradius(data&dat, int p0);
double interpolation (double x1,double x2, double y1,double y2,double x3);
double getinterval (data&dat, int p0);
// main procedure
int main(int argc, char **argv) try
	{
	string filename;
	if (argc<2)
		{
			filename = "sand.csd";
		}
	else
		{
			filename = argv[1];
		}
	data dat;
	dat.check = true;
	dat.particlereuse = false;
	dat.facereuse = false;
	dat.approximation =1.0e-9;
	int fraction=0;
	int interval=0;
	int usingparticle =0;
	int usingface = 0;
	int count=0;
	double distance1=0;
	double distance2=0;
	ifstream datain;
	//cout << "Welcome to CSD program \n";
	//cout << "Enter data file name: ";
	//cin >> dat.datafile;
	datain.open(filename.c_str());
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
	if (dat.specimentype.compare("Sphere")==0)
		{
			datain>>dat.specimen[0];
			dat.specimen[0]/=2;
			dat.volume=4./3.*M_PI*pow(dat.specimen[0],3.);
		}
	else if(dat.specimentype.compare("Cube")==0)
		{
			for (int i=0; i<3; i++)
				{
					datain>>dat.specimen[i];
					dat.specimen[i]/=2;
				}
			dat.volume=8*dat.specimen[0]*dat.specimen[1]*dat.specimen[2];
		}
	else if((dat.specimentype.compare("Cylinder")==0))
		{
			for (int i=0; i<2; i++)
				{
					datain>>dat.specimen[i];
					dat.specimen[i]/=2;
				}
			dat.volume=M_PI*pow(dat.specimen[0],2.)*dat.specimen[1];
		}
									datain.ignore(200,'\n');
	datain >> dat.approximation;	datain.ignore(200,'\n');
	datain >> dat.GSD;
	datain >> dat.numbershapes;		datain.ignore(200,'\n');
	for (int i=1; i<dat.GSD+1;i++)
		{
			datain >> fraction;
			datain >> dat.gsd[fraction];
			dat.gsd1.push_back(fraction);
			datain >> dat.densities[fraction];
			dat.densities[fraction] /=1000;																													// change to g/mm3
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
	dat.density = dat.density/1000.0;																														// transfer density from g/cm3 to g/mm3
	gsdgenerate(dat);																																		// Generate gsd
	// choose again volume
	dat.choice = false;
	while (!dat.choice)
		{
			calculateintervals(dat);																														// generate list of intervals with number of particles
			// check information after calculation
			//for (int i=0; i<dat.numberintervals; i++)
			//	{
			//		cout << i << " " << dat.intervals[i].begin<< " " << dat.intervals[i].numberparticles << " " << dat.intervals[i].diameter << "\n";
			//	}
			//cout << dat.numberparticles;
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
	listprepare(dat);																			// generate list of particles with confirmed volume
	cout << "created list of particles \n";
	//generate specimen;
	if (dat.specimentype.compare("Cube")==0)
		{
			for (int i =0; i<3; i++)
				{
					dat.specimen[i]*=dat.factor[i];
				}
			//Vec3_t dummy(dat.specimen[0],0,0);
			//dat.particles.AddPlane(-1,dummy,0.1,2*dat.specimen[0],2*dat.specimen[0],1.0,0.5*M_PI,&OrthoSys::e1);
		}
	else if (dat.specimentype.compare("Sphere")==0)
		{
			dat.specimen[0]*=dat.factor[0];
		}
	else if (dat.specimentype.compare("Cylinder")==0)
		{
			for (int i =0; i<2; i++)
					{
						dat.specimen[i]*=dat.factor[i];
					}
		}
	dat.numberunusedparticles = dat.numberparticles;
	// generate basic tetrahedron
	double r[3];
	dat.usingparticle = dat.numberparticles -1;
	dat.usinginterval = dat.numberintervals -1;
	while (dat.usinginterval == false)
		{
			dat.usinginterval -=1;
		}
	for (int i =0; i< 3; i++)
		{
			r[i]=equalradius(dat, dat.usingparticle-i);
		}
	dat.particles.Particles[dat.usingparticle]->Position(Vec3_t(0.,0.,0.));
	dat.particles.Particles[dat.usingparticle-1]->Position(Vec3_t(r[0]+r[1],0.,0.));
	double x = (pow(r[0]+r[1],2)+pow(r[0]+r[2],2)-pow(r[1]+r[2],2))/2/(r[0]+r[1]);
	dat.particles.Particles[dat.usingparticle-2]->Position(Vec3_t(x, pow(pow(r[0]+r[2],2)-pow(x,2),0.5),0));
	createface(dat, dat.usingparticle, dat.usingparticle-1, dat.usingparticle-2,dat.usingparticle-3);								//create face from 3 particles
	putparticle(dat,dat.usingparticle-3,0);																							//put first particle on first face, this face will put on twice
	createtetrahedron(dat,0,dat.usingparticle-3);
	dat.faces[0].faceuse=1;																											// mark used side
	dat.usingface=0;
	dat.numberfaces =1;
	dat.numberopenfaces =1;
	for (int i=0; i<4; i++)																							//mark used particles
		{
			useparticle(dat, dat.usingparticle - i);
		}
	dat.usingparticle -=4;
	usingparticle = dat.usingparticle;
	cout << "created base \n";
	textout(dat);
	// first loop for putting particles
	int usinginterval = dat.usinginterval;
	while ((dat.numberunusedparticles> 0)and(dat.numberopenfaces > 0))
		{
			while (dat.particlesuse[dat.usingparticle])
				{
					dat.usingparticle -=1;
					if (dat.usingparticle==0)
						{
							break;
						}
					}
			while (dat.faces[dat.usingface].faceuse == 0)
				{
					dat.usingface +=1;
					if (dat.usingface == dat.numberfaces)
						{
							break;
						}
				}
			//cout << dat.usingparticle << " "<<dat.usingface<< "\n";
			tryputparticle(dat,dat.usingparticle, dat.usingface);
			// build new tetrahedron with new faces
			if (dat.checkradius and dat.checkoverlap and dat.checkdistance)
				{

					createtetrahedron(dat, dat.usingface, dat.temproraryparticle);										// create new tetrahedron
					if (dat.particlereuse)
						{
							dat.particlereuse = false;
						}
					else
						{
							useparticle(dat, dat.temproraryparticle);													// mark particle used
			//				dat.usingparticle -=1;
						}
				}
			//			dat.intervals[usinginterval].usingparticle -=1;
			//			if (dat.intervals[usinginterval].usingparticle < dat.intervals[usinginterval].firstparticle)
			//				{
			//					dat.intervals[usinginterval].ability=false;
			//					usinginterval -=1;
			//					if (usinginterval <0)
			//						{
			//							dat.check =true;
			//						}
			//				}
				closeface(dat, dat.usingface);
			//	dat.usingface +=1;
			//	cout << dat.usingparticle << "\n";
			}
		cout << "put particles \n";
		cout << dat.numberparticles<< " "<<dat.numberunusedparticles << " " <<dat.usingparticle <<" " <<dat.numberfaces<< " "<<dat.usingface<< " " <<dat.numberopenfaces;
		count =0;
		for (int i=0; i< dat.particles.Particles.Size(); i++)
			{
				if (dat.particlesuse[i]==false)
					{
						dat.particles.Particles[i]->Tag -=100*dat.numbershapes;
						count+=1;
					}
			}
		// export results
		savedomain(dat);
		textout(dat);
		cout << count <<"real number unused \n";
		//for (size_t i =0; i<dat.particles.Particles.Size();i++)
		//	{
		//		checkoverlap(dat,i);
		//		if (!dat.checkoverlap)
		//		{
		//			cout << i<<"\n";
		//		}
		//	}

	    return 0;
	}
	MECHSYS_CATCH
// function
inline double equalradius(data&dat, int p0)
	{
		double r= dat.intervals[abs(int(dat.particles.Particles[p0]->Tag/dat.numbershapes))].equalradius[abs(dat.particles.Particles[p0]->Tag)%dat.numbershapes];
		r+=0.05*dat.intervals[getinterval(dat,p0)].radius;
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
inline void checkdistance(data&dat, int p0)
	{
		dat.checkdistance =true;
		if (dat.specimentype.compare("Cube")==0)
			{
				for (int i = 0; i < 3; i++)
					{
						if (abs(dat.particles.Particles[p0]->x(i))+equalradius(dat,p0) > dat.specimen[i])
							{
								dat.checkdistance= false;
							}
					}
			}
		else if (dat.specimentype.compare("Sphere")==0)
			{
				if (norm(dat.particles.Particles[p0]->x)+ equalradius(dat,p0)>dat.specimen[0])
					{
						dat.checkdistance=false;
					}
			}
	}
inline void checkoverlap(data&dat, int p0)
	{
		dat.checkoverlap =true;
		dat.overlappingpoint =-1;
		double overlapdistance = -dat.approximation;																												// default overlappint distance, to avoid calculation error
		double r1 =equalradius(dat,p0);
		double distance = 0.0;
		for ( int i = 0;  i < dat.numberparticles; i++)
			{
				if ((dat.particlesuse[i])and(!(i==p0)))
					{
						distance = norm(dat.particles.Particles[p0]->x-dat.particles.Particles[i]->x) -r1 - equalradius(dat,i);
						if (distance < overlapdistance)
							{
								dat.overlappingpoint = i;
								overlapdistance = distance;
								dat.checkoverlap =false;
							}
					}
			}
	}
inline void closeface(data & dat, int m0)
	{
		dat.faces[m0].faceuse =0;
		dat.numberopenfaces -=1;
	}
inline void createface(data & dat, int p0, int p1, int p2,int p3)
	{
		// sort particles by radius
		int sort[3];
		sort[0]=p0;
		sort[1]=p1;
		sort[2]=p2;
		int count;
		for (int i =0; i<2;i++)
			{
				for (int j=1; j<3; j++)
				{
					if (sort[i]< sort[j])
						{
							count= sort[i];
							sort[i]=sort[j];
							sort[j]=count;
						}
				}
			}
		// create face
		face facetemprorary;
		Array<Vec3_t> V(3);
		for (int i=0;i<3;i++)
			{
				facetemprorary.points[i]=sort[i];
				V[i]=dat.particles.Particles[sort[i]]->x;
			}
		// need to calculate constriction ?
		DEM::Face facedemtemprorary(V);															// add face in library
		facedemtemprorary.Normal(facetemprorary.normal);
		facetemprorary.normal /= norm(facetemprorary.normal);
		facetemprorary.point = p3;
		double distance = dot(dat.particles.Particles[p3]->x - dat.particles.Particles[facetemprorary.points[0]]->x, facetemprorary.normal);
		if (distance>0)
			{
				facetemprorary.faceuse = 1;
			}
		else
			{
				facetemprorary.faceuse =-1;
			}
		if (dat.particlereuse)
			{
				dat.facereuse =false;
				for (int i=0; i< dat.numberfaces; i++)
					{
						if (facetemprorary.points[0]== dat.faces[i].points[0])
							{
								if ((facetemprorary.points[1]== dat.faces[i].points[1])and(facetemprorary.points[1]== dat.faces[i].points[1]))
									{
										if (facetemprorary.faceuse*dat.faces[i].faceuse <0)
											{
												closeface(dat,i);
											}
										else
											{
												dat.check =true;
											}
										dat.facereuse=true;
										dat.facereusenumber = i;
									}
							}
					}
			}
		if (dat.facereuse)
			{
			dat.facereuse = false;
			}
		else
		{
			dat.faces.push_back(facetemprorary);												// add user-defined face
			dat.numberfaces +=1;
			dat.numberopenfaces +=1;
		}
	}
inline void createtetrahedron(data & dat, int m0, int p0)										// no need to create tetrahedron at the first step because the tetrahedron net must be generated before the second step
	{
		tetrahedron temprorarytetrahedron;
		temprorarytetrahedron.faces[0]= m0;
		createface(dat, dat.faces[m0].points[0], dat.faces[m0].points[1], p0, dat.faces[m0].points[2]);
		if (dat.facereuse)
			{
				temprorarytetrahedron.faces[1]= dat.facereusenumber;
				dat.facereuse=false;
			}
		else
			{
				temprorarytetrahedron.faces[1]= dat.numberfaces-1;
			}
		createface(dat, dat.faces[m0].points[0], dat.faces[m0].points[2], p0, dat.faces[m0].points[1]);
		if (dat.facereuse)
			{
				temprorarytetrahedron.faces[2]= dat.facereusenumber;
				dat.facereuse=false;
			}
		else
			{
				temprorarytetrahedron.faces[2]= dat.numberfaces-1;
			}
		createface(dat, dat.faces[m0].points[1], dat.faces[m0].points[2], p0, dat.faces[m0].points[0]);
		if (dat.facereuse)
			{
				temprorarytetrahedron.faces[3]= dat.facereusenumber;
				dat.facereuse=false;
			}
		else
			{
				temprorarytetrahedron.faces[3]= dat.numberfaces-1;
			}
		dat.tetrahedrons.push_back(temprorarytetrahedron);
	}

inline void establishlocalsystem(data&dat, int m0)
	{
		dat.localsystem[0]= dat.particles.Particles[dat.faces[m0].points[1]]->x - dat.particles.Particles[dat.faces[m0].points[0]]->x;
		dat.localsystem[0]= dat.localsystem[0]/norm(dat.localsystem[0]);							// e1 vector
		dat.localsystem[2]= dat.faces[m0].normal;													// e3 vector
		dat.localsystem[1]= -cross(dat.localsystem[0],dat.localsystem[2]);
		dat.localsystem[1]= dat.localsystem[1]/norm(dat.localsystem[1]);							// may be not in need
		dat.localroot = dat.particles.Particles[dat.faces[m0].points[0]]->x;
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
inline void listprepare(data & dat)																												// add particles to domain
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
inline void moveon(data & dat, int p0,int m0)
	{
		Vec3_t distance = dat.particles.Particles[p0]->x - dat.particles.Particles[dat.overlappingpoint]->x;
		Vec3_t projection = dot(distance, dat.faces[m0].normal);
		double move = sqrt(pow(equalradius(dat,p0) + equalradius(dat, dat.overlappingpoint),2.0) - pow(norm(distance),2.0)+pow(norm(projection),2.0))-norm(projection);
		dat.finalposition = dat.particles.Particles[p0]->x - dat.faces[m0].faceuse*move*dat.faces[m0].normal;
		dat.particles.Particles[p0]->Position(dat.finalposition);
	}
inline void putparticle(data & dat, int p0, int m0)
	{
		establishlocalsystem(dat,m0);															// calculate local coordinates system
		dat.checkradius =false;
		double r[4];
		for (int i =0;i<3;i++)
			{
				r[i]= equalradius(dat, dat.faces[m0].points[i]);
			}
		r[3] = equalradius(dat,p0);
		double x2 = norm(dat.particles.Particles[dat.faces[m0].points[1]]->x -dat.particles.Particles[dat.faces[m0].points[0]]->x);
		double a1 = (r[0] -r[1])/x2;
		double b1 = (pow(x2,2.)+pow(r[0],2.)-pow(r[1],2.))/2./x2;
		double x3 = dot(dat.particles.Particles[dat.faces[m0].points[2]]->x - dat.particles.Particles[dat.faces[m0].points[0]]->x,dat.localsystem[0]);
		double y3 = dot(dat.particles.Particles[dat.faces[m0].points[2]]->x - dat.particles.Particles[dat.faces[m0].points[0]]->x, dat.localsystem[1]);
		double a2 = (r[0]-r[2]-a1*x3)/y3;
		double b2 = (pow(x3,2.)+pow(y3,2.)+pow(r[0],2.)-pow(r[2],2.)-2.*b1*x3)/2./y3;
		double x4 = a1*r[3]+b1;
		double y4 = a2*r[3]+b2;
		double z4 = pow(r[0]+r[3],2)-pow(a1*r[3]+b1,2)-pow(a2*r[3]+b2,2);
		if (z4 > 0)
			{
				z4=-dat.faces[m0].faceuse*pow(z4,0.5);
				Vec3_t position(x4,y4,z4);
				Vec3_t finalposition =  dat.particles.Particles[dat.faces[m0].points[0]]->x;
				for (int i=0; i<3; i++)
					{
						for (int j=0; j<3; j++)
							{
								finalposition(i)+= position(j)*dat.localsystem[j](i);
							}
					}
				dat.particles.Particles[p0]->Position(finalposition);
				dat.checkradius =true;
			}
	}
inline void savedomain(data&dat)
	{
		dat.particles.WriteXDMF(dat.outputdomain.c_str());// export to draw visual results
		dat.particles.Save(dat.outputdomain.c_str());
		cout << "saved domain: " << dat.outputdomain <<"\n";
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
inline void tryputparticle(data&dat,int p0, int m0)
	{
		dat.temproraryparticle = p0;
		int usingface = m0;
		double distancepoint =0;
		double distanceoverlap =0;
		putparticle(dat,dat.temproraryparticle,usingface);
		checkoverlap(dat, dat.temproraryparticle);
		checkdistance(dat, dat.temproraryparticle);
		if ((!dat.checkradius)or(!dat.checkdistance)or(!dat.checkoverlap))
			{
				int interval = getinterval(dat,dat.temproraryparticle);
				distancepoint = dot(dat.particles.Particles[dat.faces[usingface].point]->x - dat.particles.Particles[dat.faces[usingface].points[0]]->x, dat.faces[usingface].normal);
				if (equalradius(dat,dat.temproraryparticle) > dat.faces[usingface].constrictionsize)										// if particle size is bigger than constriction size
					{
						while ((!dat.checkdistance)or(!dat.checkradius)or(!dat.checkoverlap))												// try to reduce size of particles
							{
								interval -=1;
								if ((interval <0)or (dat.intervals[interval].radius < dat.faces[usingface].constrictionsize))
									{
										if (interval < dat.numberintervals-1)
											{
												interval+=1;
											}
										break;
									}
								if (dat.intervals[interval].ability)
									{
										dat.temproraryparticle = dat.intervals[interval].usingparticle;
										putparticle(dat, dat.temproraryparticle, usingface);
										checkdistance (dat, dat.temproraryparticle);
										checkoverlap(dat, dat.temproraryparticle);
										if ((dat.checkradius) and (dat.checkdistance) and (!dat.checkoverlap))
											{
												distanceoverlap = dot(dat.particles.Particles[dat.overlappingpoint]->x - dat.particles.Particles[dat.faces[usingface].points[0]]->x, dat.faces[usingface].normal);
												while ((distancepoint*distanceoverlap>0)and (!dat.checkoverlap)and(dat.checkdistance))							// if overlap with particle from other side of face then move particle
													{
														moveon(dat,dat.temproraryparticle,usingface);
														checkoverlap(dat, dat.temproraryparticle);
														checkdistance(dat,dat.temproraryparticle);
														if (!dat.checkoverlap)
															{
																distanceoverlap = dot(dat.particles.Particles[dat.overlappingpoint]->x - dat.particles.Particles[dat.faces[usingface].points[0]]->x, dat.faces[usingface].normal);
															}
													}
											}
									}
							}
					}
				if ((!dat.checkradius)or(!dat.checkdistance)or(!dat.checkoverlap))													// if particle size is smaller than constriction size
					{
						while ((!dat.checkdistance)or(!dat.checkradius)or(!dat.checkoverlap))										// try to reduce size of particles
							{
								interval -=1;
								if (interval <0)
									{
										break;
									}
								if (dat.intervals[interval].ability)
									{
										dat.temproraryparticle = dat.intervals[interval].usingparticle;
										dat.particles.Particles[dat.temproraryparticle]->Position(dat.faces[usingface].constrictioncentre-dat.approximation*dat.faces[usingface].normal*dat.faces[usingface].faceuse);						// put particle at the centre of constriction
										checkdistance (dat, dat.temproraryparticle);
										checkoverlap(dat, dat.temproraryparticle);
										if ((dat.checkradius) and (dat.checkdistance) and (!dat.checkoverlap))
											{
												distanceoverlap = dot(dat.particles.Particles[dat.overlappingpoint]->x - dat.particles.Particles[dat.faces[usingface].points[0]]->x, dat.faces[usingface].normal);
												while ((distancepoint*distanceoverlap>0)and (!dat.checkoverlap)and(dat.checkdistance))							// if overlap with particle from other side of face then move particle
													{
														moveon(dat,dat.temproraryparticle,usingface);
														checkoverlap(dat, dat.temproraryparticle);
														checkdistance(dat,dat.temproraryparticle);
														if (!dat.checkoverlap)
															{
																distanceoverlap = dot(dat.particles.Particles[dat.overlappingpoint]->x - dat.particles.Particles[dat.faces[usingface].points[0]]->x, dat.faces[usingface].normal);
															}
													}
											}
									}
							}
					}
			}
	}
inline void useparticle(data & dat, int p0)
	{
		if (dat.particlesuse[p0])
			{
				cout << "has been used before \n";
			}
		dat.particlesuse[p0]=true;
		dat.numberunusedparticles -=1;
		int interval = abs(int(dat.particles.Particles[p0]->Tag/dat.numbershapes));
		dat.intervals[interval].usingparticle -=1;
		if (dat.intervals[interval].usingparticle < dat.intervals[interval].firstparticle )
			{
				dat.intervals[interval].ability =false;
				if (interval == dat.usinginterval)
					{
						dat.usinginterval -=1;
					}
			}

	}
