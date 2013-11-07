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
	int faceuse;
	int points[3];
	int point;
	int overlaps[3];
	double constrictionsize;
	double smallsizes[3];
	Vec3_t centre;
	Vec3_t constrictioncentre;
	Vec3_t normal;
	Vec3_t smallcentres[3];
	Vec3_t systems[3];
	bool reduce;
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
	int overlappingpoint;
	int usingface;
	int usinginterval;
	int usingparticle;
	int temproraryparticle;
	string datafile;
	string outputdomain;
	string soilname;
	string line;
	string specimentype;
	string fillingtype;
	string particleshapetype;
	double specimen[3];
	double factor[3];
	double gsd[101];
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
void checkface(data&dat,int m0=1001);
void checkparticle(data&dat,int p0);
void closeface(data&dat, int m0);
void convertcoordinates(data & dat, Vec3_t & position, bool way);
void createface(data & dat, int p0, int p1, int p2,int faceuse);
void createtetrahedron(data & dat, int m0, int p0);
void establishlocalsystem(data&dat, int m0);
void gsdgenerate(data & dat);
void listprepare(data & dat);
void moveon(data & dat, int p0, int m0);
void reduceconstriction( data&dat, int m0, int overlap);
void planeconstriction(data&dat, int m0, int p0, int p1, int overlap, int smallconstriction);
void putparticle(data & dat, int p0, int m0);
void savedomain(data&dat);
void textout (data& dat);
void textoutconstrictionsize(data&dat);
void tryputparticle(data&dat, int p0,int m0);
void useparticle(data & dat, int p0);
double interpolation (double x1,double x2, double y1,double y2,double x3);
double constrictionsize(data&dat, int p0, int p1, int p2);
int checksphereoverlap(data&dat, double radius, Vec3_t centre);
// main procedure
int main(int argc, char**argv)try
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
	int fraction=0;
	int count =0;
	ifstream datain;
	cout << "Reading data \n";
	datain.open(filename.c_str());
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
	datain >> dat.GSD;				datain.ignore(200,'\n');
	for (int i=1; i<dat.GSD+1;i++)
		{
			datain >> fraction;
			datain >> dat.gsd[fraction];
			dat.gsd1.push_back(fraction);
			datain.ignore (200,'\n');
		}
	for (int i=0; i<3; i++)
		{
			datain >> dat.factor[i];
		}
									datain.ignore (200,'\n');
	datain >> dat.outputdomain;		datain.ignore (200,'\n');
	datain.close();
	dat.density = dat.density/1000.0;																			// transfer density from g/cm3 to g/mm3
	gsdgenerate(dat);																																		// Generate gsd
	// choose again volume
	dat.choice = false;
	cout << "Calculating intervals \n";
	while (!dat.choice)
		{
			calculateintervals(dat);
			// generate list of intervals with number of particles
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
	cout << "Preparing list \n";
	listprepare(dat);																			// generate list of particles with confirmed volume
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
	textout (dat);
	// generate basic tetrahedron
	cout << "Generating basic tetrahedron \n";
	double r[3];
	dat.usingparticle = dat.numberparticles -1;
	dat.usinginterval = dat.numberintervals -1;
	for (int i =0; i< 3; i++)
		{
			r[i]=dat.particles.Particles[dat.usingparticle-i]->Props.R;
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
	// first loop for putting particles
	cout << "Sequential packing \n";
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
					dat.usingface+=1;
					if (dat.usingface == dat.numberfaces-1)
						{
							break;
						}
				}
			tryputparticle(dat,dat.usingparticle, dat.usingface);
			// build new tetrahedron with new faces
			if (dat.checkradius and dat.checkoverlap and dat.checkdistance)
				{
					closeface(dat, dat.usingface);
					createtetrahedron(dat, dat.usingface, dat.temproraryparticle);										// create new tetrahedron
					if (dat.particlereuse)
						{
							dat.particlereuse = false;
							//cout << "hi \n";
						}
					else
						{
							useparticle(dat, dat.temproraryparticle);													// mark particle used
						}
					//dat.intervals[usinginterval].usingparticle -=1;
					//if (dat.intervals[usinginterval].usingparticle < dat.intervals[usinginterval].firstparticle)
					//	{
					//		dat.intervals[usinginterval].ability=false;
					//		usinginterval -=1;
					//	}
				}
			else
				{
					closeface(dat, dat.usingface);
				}
		}
	for (int i=0; i< dat.particles.Particles.Size(); i++)
		{
			if (dat.particlesuse[i]==false)
			{
				dat.particles.Particles[i]->Tag -=100;
			}

		}
	cout << dat.numberparticles<< " "<<dat.numberunusedparticles << " " <<dat.usingparticle <<" " <<dat.numberfaces<< " "<<dat.usingface<< " " <<dat.numberopenfaces<<"\n";
	dat.check=false;
	for (int i=0;i<dat.numberparticles;i++)
		{
			if (dat.particlesuse[i]==false)
				{
					dat.particles.Particles[i]->Tag = -1000;
					dat.check =true;
				}
		}
	if (dat.check)
	{
		Array <int> delpar;
		delpar.Push(-1000);
		cout <<delpar[0] << "\n";
		cout << "deleting \n";
		dat.particles.DelParticles(delpar);
	}
	// export results
	savedomain(dat);
	textout (dat);
	// constriction calculation
	dat.count=0;
	cout << "Recalculating constriction \n";
	int overlap;
	for (int i=0; i <dat.faces.size(); i++)
	{
//		if (dat.faces[i].constrictionsize==0.)
//		{
//			continue;
//		}
		if (dat.faces[i].constrictionsize<=0.)
		{
			cout << "error \n";
		}
		overlap =checksphereoverlap(dat, dat.faces[i].constrictionsize, dat.faces[i].constrictioncentre);
		if (overlap >-1)
		{
			bool check=false;
			for (int j=0; j<3;j++)
			{
				if (overlap ==dat.faces[i].points[j])
				{
//					cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n";
//					checkface(dat,i);
//					checkparticle(dat,overlap);
//					double distance = norm (dat.particles.Particles[overlap]->x-dat.faces[i].constrictioncentre)-dat.particles.Particles[overlap]->Props.R-dat.faces[i].constrictionsize;
//					cout << distance << "\n";
					check=true;
				}
			}
			if (check)
			{
				continue;
			}
			dat.faces[i].reduce=true;
			reduceconstriction(dat,i,overlap);
		}
		else
		{
			dat.faces[i].reduce=false;
		}
	}
	textoutconstrictionsize(dat);
}
MECHSYS_CATCH
// function
inline double interpolation (double x1,double x2, double y1,double y2,double x3)
	{
		double y3;
		y3= y1 + (x3-x1)*(y2-y1)/(x2-x1);
		return (y3);
	}
// procedure
inline double constrictionsize (data&dat, int p0, int p1, int p2)																	// this constriction size is not real, because it may overlap with particles
	{
		Vec3_t x = dat.particles.Particles[p1]->x-dat.particles.Particles[p0]->x;
		double x2 = norm(x);																										// get coordinates of particles in local coordinate system
		x /=x2;
		double x3 = dot(dat.particles.Particles[p2]->x - dat.particles.Particles[p0]->x, x);
		Vec3_t y = (dat.particles.Particles[p2]->x-dat.particles.Particles[p0]->x-x3*x);
		double y3 = norm(y);
		y=y/y3;
		double a1 = (dat.particles.Particles[p0]->Props.R-dat.particles.Particles[p1]->Props.R)/x2;
		double a2 = (dat.particles.Particles[p0]->Props.R -dat.particles.Particles[p2]->Props.R-a1*x3)/y3;
		double b1 = (x2*x2+pow(dat.particles.Particles[p0]->Props.R,2.0)-pow(dat.particles.Particles[p1]->Props.R,2.0))/2.0/x2;
		double b2 = (pow(x3,2.0)+pow(y3,2.0)+pow(dat.particles.Particles[p0]->Props.R,2.0)-pow(dat.particles.Particles[p2]->Props.R,2.0)-2*b1*x3)/2.0/y3;
		double a = pow(a1,2.0)+pow(a2,2.0)-1;
		double b = a1*b1+a2*b2-dat.particles.Particles[p0]->Props.R;
		double c = pow(b1, 2.0) +pow(b2, 2.0) - pow(dat.particles.Particles[p0]->Props.R, 2.0);
		if (b*b-a*c>0.)
		{
			double constriction = (-b-pow(pow(b, 2.0)-a*c,0.5))/a;
			if (constriction >dat.approximation)
			{
				dat.position = (a1*constriction+b1)*x+(a2*constriction+b2)*y+dat.particles.Particles[p0]->x;
				return (constriction);
			}
			else
			{
				cout << "constriction problems \n";
				return (0.0);
			}
		}
		else
		{
			cout << "no constriction \n";
			cout << x<< " "<< y<< " \n";
			cout << x2 << " " << x3<< " "<<y3 <<"\n";
			cout << a1 << " " << a2 <<" "<< b1<< " "<< b2<<" \n";
			cout << a << " " << b << " " << c<< "\n";
			checkparticle(dat,p0);
			checkparticle(dat,p1);
			checkparticle(dat,p2);
			return(0.0);
		}
	}
inline int checksphereoverlap(data & dat, double radius, Vec3_t constrictioncentre)
	{
		bool overlapcheck = true;
		int overlap;
		double overlapdistance = -dat.approximation;
		double distance;
		for (int i =0; i < dat.particles.Particles.Size();i++)
		{
			distance = norm(dat.particles.Particles[i]->x-constrictioncentre) -(dat.particles.Particles[i]->Props.R+radius);
			if (distance < overlapdistance)
			{
				overlapcheck =false;
				overlapdistance = distance;
				overlap =i;
				dat.count +=1;
			}
//		if (abs(dot(dat.faces[m0].constrictioncentre-dat.particles.Particles[dat.faces[m0].points[0]]->x,dat.faces[m0].normal))>0.0001)
//		{
//			cout << "constriction out of plane \n";
//		}
		}
		if (overlapcheck)
		{
			overlap =-1;
		}
		return(overlap);
	}
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
		double mass = dat.density * dat.volume * (1-dat.porosity);																							// if density of soil is constant
		double module =0;																																// abundant mass
		dat.numberparticles =0;
		for (int i =dat.numberintervals-1; i> -1; i--)
			{
				dat.intervals[i].end = dat.intervals[i+1].begin;																						// not completely delete before
				dat.intervals[i].begindiameter = dat.gsd[dat.intervals[i].begin];
				dat.intervals[i].enddiameter = dat.gsd[dat.intervals[i].end];
				dat.intervals[i].density= dat.density; 																									// if density of soil is constant
				dat.intervals[i].mass = mass*(dat.intervals[i].end-dat.intervals[i].begin)/100.0 + module;
				dat.intervals[i].diameter = pow(10.0, (log10(dat.intervals[i].begindiameter)+log10(dat.intervals[i].enddiameter))/2); 					// generate mean diameter of interval
				dat.intervals[i].radius = dat.intervals[i].diameter/2.0;
				double particlemass = dat.intervals[i].density * pow(dat.intervals[i].diameter,3.0)/3.0*acos(0.0);									// mass of one particle must be transfer to g/cm3
				dat.intervals[i].numberparticles = dat.intervals[i].mass/particlemass;														 			// calculate number of particles, pi = 2*acos(0)
				module = dat.intervals[i].mass-dat.intervals[i].numberparticles*particlemass;															// abundant mass is passed to next interval
				dat.intervals[i].mass = dat.intervals[i].numberparticles*particlemass;																	// real mass of interval after calculation
				dat.numberparticles += dat.intervals[i].numberparticles;
			}
	}
inline void checkdistance(data&dat, int p0)
	{
		dat.checkdistance =true;
		if (dat.specimentype.compare("Cube")==0)
			{
				for (int i = 0; i < 3; i++)
					{
						if (abs(dat.particles.Particles[p0]->x(i))+dat.particles.Particles[p0]->Props.R > dat.specimen[i])
							{
								dat.checkdistance= false;
							}
					}
			}
		else if (dat.specimentype.compare("Sphere")==0)
			{
				if (norm(dat.particles.Particles[p0]->x)+dat.particles.Particles[p0]->Props.R>dat.specimen[0])
					{
						dat.checkdistance=false;
					}
			}
	}

inline void checkface(data&dat,int m0)
	{
		for (int i=0; i<3;i++)
		{
			cout << dat.faces[m0].points[i]<< " "<< dat.particles.Particles[dat.faces[m0].points[i]]->x <<" " <<dat.particles.Particles[dat.faces[m0].points[i]]->Props.R << "\n";
		}
		cout << dat.faces[m0].constrictioncentre<< " "<<dat.faces[m0].constrictionsize << "\n";
	}
inline void checkoverlap(data&dat, int p0)
	{
		dat.checkoverlap =true;
		dat.overlappingpoint =-1;
		double overlapdistance = -dat.approximation;																												// default overlappint distance, to avoid calculation error
		for ( int i = 0;  i < dat.numberparticles; i++)
			{
				if ((dat.particlesuse[i])and(!(i==p0)))
					{
						double distance = norm(dat.particles.Particles[p0]->x-dat.particles.Particles[i]->x) -dat.particles.Particles[p0]->Props.R -dat.particles.Particles[i]->Props.R;
						if (distance < overlapdistance)
							{
								dat.overlappingpoint = i;
								overlapdistance = distance;
								dat.checkoverlap =false;
							}
					}
			}
	}
inline void checkparticle(data&dat, int p0)
{
	cout << p0<< " "<< dat.particles.Particles[p0]->x<< " "<< dat.particles.Particles[p0]->Props.R <<" \n";
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
		if (dat.particlereuse)																																	// check whether face is reused or not
					{
						dat.facereuse =false;
						for (int i=0; i< dat.numberfaces; i++)
							{
								if (facetemprorary.points[0]== dat.faces[i].points[0])
									{
										if ((facetemprorary.points[1]== dat.faces[i].points[1])and(facetemprorary.points[2]== dat.faces[i].points[2]))
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
				facetemprorary.constrictionsize = constrictionsize(dat, facetemprorary.points[0],facetemprorary.points[1],facetemprorary.points[2]);		// calculate constriction size
				facetemprorary.constrictioncentre = dat.position;
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
						}
				}
	}
inline void listprepare(data & dat)																								// add particles to domain
	{
		int count =0;
		Vec3_t x(0.,0.,0.);
		for (int i=0;i<dat.numberintervals;i++)
			{
				dat.intervals[i].firstparticle = count;
				for (int j=0; j<dat.intervals[i].numberparticles; j++)
					{
						dat.particles.AddSphere(i, x ,dat.intervals[i].radius, dat.intervals[i].density); 				//e - coordinates has been declared before
						dat.particlesuse.push_back(false);
						count +=1;
					}
				dat.intervals[i].lastparticle =count -1;
				dat.intervals[i].usingparticle =count-1;
				dat.intervals[i].ability =true;
			}
	}
inline void moveon(data & dat, int p0, int m0)
	{
		Vec3_t distance = dat.particles.Particles[p0]->x - dat.particles.Particles[dat.overlappingpoint]->x;
		Vec3_t projection = dot(distance, dat.faces[m0].normal);
		double move = sqrt(pow(dat.particles.Particles[p0]->Props.R + dat.particles.Particles[dat.overlappingpoint]->Props.R,2.0) - pow(norm(distance),2.0)+pow(norm(projection),2.0))-norm(projection);
		dat.finalposition = dat.particles.Particles[p0]->x - dat.faces[m0].faceuse*move*dat.faces[m0].normal;
		dat.particles.Particles[p0]->Position(dat.finalposition);
	}
inline void reduceconstriction(data&dat, int m0, int overlap)
	{
		establishlocalsystem(dat,m0);
		int checkID=overlap;
		int points[10];
		for (int i=0; i<3; i++)
		{
			int count =0;
			while (checkID>-1)
			{
				points[count]=checkID;
				//cout << checkID << "\n";
				planeconstriction(dat, m0, dat.faces[m0].points[i], dat.faces[m0].points[(i+1)%3], checkID, i);
				checkID=checksphereoverlap(dat,dat.faces[m0].smallsizes[i],dat.faces[m0].smallcentres[i]);
				if ((checkID==dat.faces[m0].points[i])or(checkID==dat.faces[m0].points[(i+1)%3]))
				{
					cout << "overlap with particle on face \n";
					checkface(dat,m0);
					checkparticle(dat,checkID);
					cout << dat.faces[m0].smallcentres[i]<< " "<<dat.faces[m0].smallsizes[i]<<" \n";
				}
				if (checkID==dat.faces[m0].points[(i+2)%3])
				{
					double distance = norm (dat.particles.Particles[checkID]->x-dat.faces[m0].smallcentres[i])-dat.particles.Particles[checkID]->Props.R-dat.faces[m0].smallsizes[i];
					if (distance <dat.particles.Particles[checkID]->Props.R/20.)
					{
						checkID=-1;
					}
				}
				count +=1;
				if (count >9)
				{
//					cout <<" long circles \n";
//					for (int j=0; j <10; j++)
//					{
//						cout << points[j]<< " ";
//					}
//					cout <<"\n";
					checkID=-1;
					dat.faces[m0].smallsizes[i]=0.;
					break;
				}
			}
		}
	}
inline void planeconstriction(data&dat, int m0, int p0, int p1, int overlap, int smallconstriction)
{
	Vec3_t system[3];
	int p[3];
	p[0]=p0;
	p[1]=p1;
	p[2]=overlap;
	int swap;
	if (dat.particles.Particles[p[0]]->Props.R < dat.particles.Particles[p[1]]->Props.R)
	{
		swap = p[0];
		p[0]=p[1];
		p[1]=swap;
	}
	system[2]=dat.localsystem[2];
	system[0]= (dat.particles.Particles[p[1]]->x-dat.particles.Particles[p[0]]->x);
	double x2 = norm(system[0]);
	system[0]/=x2;
	system[1]=-cross(system[0],system[2]);
	system[1]/=norm(system[1]);
	double x3[3];
	double r[3];
	for (int i=0;i<3;i++)
	{
		x3[i]=dot(dat.particles.Particles[p[2]]->x-dat.particles.Particles[p[0]]->x,system[i]);
		r[i]=dat.particles.Particles[p[i]]->Props.R;
	}
	double a1 = (r[0]-r[1])/x2;
	double a2 = (r[0]*r[0]+x2*x2-r[1]*r[1])/2.0/x2;
	if (x2==0)
	{
		cout << "fault \n";
	}
	double b1 = (r[0]-r[2]-a1*x3[0])/x3[1];
	double b2 = (r[0]*r[0]-r[2]*r[2]+x3[0]*x3[0]+x3[1]*x3[1]+x3[2]*x3[2]-2*a2*x3[0])/2.0/x3[1];
	if (x3[1]==0)
	{
		cout << "fault2 \n";
		cout << p[0]<<" "<<p[2];
	}
	double a = a1*a1+b1*b1-1.0;
	double b = a1*a2+b1*b2-r[0];
	double c = a2*a2+b2*b2-r[0]*r[0];
//	double x[3];
//	double y[3];
//	double z[3];
//	double r[3];
//	for (int i=0; i<3; i++)
//	{
//		x[i]=dot(dat.particles.Particles[p[i]]->x-dat.particles.Particles[dat.faces[m0].points[0]]->x,dat.localsystem[0]);
//		y[i]=dot(dat.particles.Particles[p[i]]->x-dat.particles.Particles[dat.faces[m0].points[0]]->x,dat.localsystem[1]);
//		z[i]=dot(dat.particles.Particles[p[i]]->x-dat.particles.Particles[dat.faces[m0].points[0]]->x,dat.localsystem[2]);
//		r[i]=dat.particles.Particles[p[i]]->Props.R;
//	}
//	double num;
//	num = pow(x[0],2.)+pow(y[0],2.)+pow(z[0],2.);
//	double a1 = (pow(r[1],2.)-pow(r[0],2.)+num-(pow(x[1],2.)+pow(y[1],2.)+pow(z[1],2.)))/2;
//	double a2 = (pow(r[2],2.)-pow(r[0],2.)+num-(pow(x[2],2.)+pow(y[2],2.)+pow(z[2],2.)))/2;
//	num = (y[0]-y[2])*(x[0]-x[1])-(y[0]-y[1])*(x[0]-x[2]);
//	if (num==0)
//	{
//		cout << x[0]<< " "<< x[1]<< " "<< x[2]<< " " <<y[0]<< " " << y[1] << " " <<y[2]<<" fault \n";
//	}
//	double b1 = ((r[2]-r[0])*(x[0]-x[1])-(r[1]-r[0])*(x[0]-x[2]))/num;
//	double b2 = (a2*(x[0]-x[1])-a1*(x[0]-x[2]))/num;
//	double c1 = (r[1]-r[0]-(y[0]-y[1])*b1)/(x[0]-x[1]);
//	double c2 = (a1-(y[0]-y[1])*b2)/(x[0]-x[1]);
//	if (x[0]==x[1])
//	{
//		cout << x[0]<< " "<< x[1]<< " "<< x[2]<< " " <<y[0]<< " " << y[1] << " " <<y[2]<<" fault2 \n";
//	}
//	double a = pow(c1,2.)+pow(b1,2.)-1;
//	double b = c1*(c2-x[0])+b1*(b2-y[0])-r[0];
//	double c = pow(c2-x[0],2.)+pow(b2-y[0],2.)-r[0]*r[0]+z[0]*z[0];
	if (b*b>a*c)
	{
		double constriction = (-b-pow(b*b-a*c,0.5))/a;
		if ((constriction >dat.approximation)and(x3[1]*(b1*constriction+b2)>0))
		{
			dat.faces[m0].smallsizes[smallconstriction]=constriction;
			dat.faces[m0].smallcentres[smallconstriction]=(a1*constriction+a2)*system[0]+(b1*constriction+b2)*system[1]+dat.particles.Particles[p[0]]->x;
			//dat.faces[m0].smallcentres[smallconstriction]=dat.particles.Particles[dat.faces[m0].points[0]]->x+(c1*constriction+c2)*dat.localsystem[0]+(b1*constriction+b2)*dat.localsystem[1];
//			for (int i=0;i<3;i++)
//			{
//				cout << p[i]<< " "<< dat.particles.Particles[p[i]]->x<< " "<<dat.particles.Particles[p[i]]->Props.R<< "\n";
//			}
//			cout << dat.faces[m0].smallcentres[smallconstriction]<< " "<<dat.faces[m0].smallsizes[smallconstriction]<<"\n";
		}
		else
		{
			constriction = (-b+pow(b*b-a*c,0.5))/a;
			if (constriction >dat.approximation)
			{
				dat.faces[m0].smallsizes[smallconstriction]=constriction;
				dat.faces[m0].smallcentres[smallconstriction]=(a1*constriction+a2)*system[0]+(b1*constriction+b2)*system[1]+dat.particles.Particles[p[0]]->x;
				//dat.faces[m0].smallcentres[smallconstriction]=dat.particles.Particles[dat.faces[m0].points[0]]->x+(c1*constriction+c2)*dat.localsystem[0]+(b1*constriction+b2)*dat.localsystem[1];
				//cout << " swap sign constriction \n";
			}
			else
			{
//				cout << p[0] <<" "<< p[1] <<" "<<p[2]<<" small constriction problem \n";
//				checkface(dat,m0);
//				checkparticle(dat,overlap);
				dat.faces[m0].smallsizes[smallconstriction]=0.;
			}
		}
	}
	else
	{
		cout << "no small constriction \n";
	}
}
inline void putparticle(data & dat, int p0, int m0)
	{
		establishlocalsystem(dat,m0);															// calculate local coordinates system
		dat.checkradius =false;
		double r[4];
		for (int i =0;i<3;i++)
			{
				r[i] = dat.particles.Particles[dat.faces[m0].points[i]]->Props.R;
			}
		r[3]=dat.particles.Particles[p0]->Props.R;
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
inline void textoutconstrictionsize(data&dat)
	{
		vector <double> constrictionsize;							// get all constrictions
		for (size_t i=0; i <dat.faces.size();i++)
		{
			if (!dat.faces[i].reduce)
			{
				if (dat.faces[i].constrictionsize>dat.approximation)
				{
					constrictionsize.push_back(dat.faces[i].constrictionsize);
				}
			}
			else
			{
				for (int j=0; j<3;j++)
				{
					if (dat.faces[i].smallsizes[j]>dat.approximation)
					{
						constrictionsize.push_back(dat.faces[i].smallsizes[j]);
					}
				}
			}
		}
		double max = constrictionsize[0];							// define log base
		double min = constrictionsize[0];
		for (size_t i=0;i<constrictionsize.size();i++)
		{
			if (max<constrictionsize[i])
			{
				max=constrictionsize[i];
			}
			if (min>constrictionsize[i])
			{
				min = constrictionsize[i];
			}
		}
		cout << min << " " << max<<"\n";
		int powermin = floor(log10(min));
		int powermax = floor (log10(max));
		int count = floor(min/pow(10.,powermin));
		double power=powermin;
		vector <double> constriction;
		vector <double> number;
		vector <double> area;
		vector <double> mass;
		double siz=count*pow(10,power);
		double zer =0.;
		cout <<zer <<"\n";
		while (siz<max)
		{
			siz=count*pow(10.,power);
			constriction.push_back(siz);
			cout <<zer <<"\n";
			number.push_back(zer);
			area.push_back(zer);
			mass.push_back(zer);
			count =(count+1)%10;
			if (count==0)
			{
				power+=1;
				count=1;
			}
		}
		constriction.push_back(count*pow(10.,power));
		number.push_back(zer);
		area.push_back(zer);
		mass.push_back(zer);
		cout << "accumulation \n";
		for (int i=1;i<constriction.size();i++)									// accumulate fraction
		{
			for (size_t j=0; j<constrictionsize.size();j++)
			{
				if (constrictionsize[j]<constriction[i])
				{
					number[i]+=1;
					area[i]+=acos(-1.)*pow(constrictionsize[j],2.);
					mass[i]+=4/3*acos(-1.)*pow(constrictionsize[j],3.);
				}
			}
			if (i<constriction.size()-1)
			{
				number[i+1]=number[i];
				area[i+1]=area[i];
				mass[i+1]=mass[i];
			}
		}
		for (int i=0;i<constriction.size();i++)									// normalise
		{
			number[i]=number[i]*100/number[constriction.size()-1];
			area[i]=area[i]*100/area[constriction.size()-1];
			mass[i]=mass[i]*100/mass[constriction.size()-1];
		}
		ofstream dataout;																							// export to text file
		dataout.open("csd.out");
		dataout << "size number area mass \n";
		for (int i=0; i<constriction.size();i++)
		{
			dataout<< constriction[i]<< " "<<number[i]<<" "<< area[i]<< " "<< mass[i]<<"\n";
		}
		dataout.close();
	}
inline void tryputparticle(data&dat,int p0, int m0)
	{
		dat.temproraryparticle = p0;
		int usingface = m0;
		double distancepoint =0;
		double distanceoverlap =0;
		putparticle(dat, dat.temproraryparticle,usingface);
		checkoverlap(dat, dat.temproraryparticle);
		checkdistance(dat, dat.temproraryparticle);
		if ((!dat.checkradius)or(!dat.checkdistance)or(!dat.checkoverlap))
			{
				int interval = dat.particles.Particles[dat.temproraryparticle]->Tag;
				distancepoint = dot(dat.particles.Particles[dat.faces[usingface].point]->x - dat.particles.Particles[dat.faces[usingface].points[0]]->x, dat.faces[usingface].normal);
				if (dat.particles.Particles[dat.temproraryparticle]->Props.R > dat.faces[usingface].constrictionsize)			// if particle size is bigger than constriction size
					{
						while ((!dat.checkdistance)or(!dat.checkradius)or(!dat.checkoverlap))										// try to reduce size of particles
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
						while ((!dat.checkdistance)or(!dat.checkoverlap)or(!dat.checkradius))										// try to reduce size of particles
							{
								interval -=1;
								if (interval <0)
									{
										break;
									}
								if (dat.intervals[interval].ability)
									{
										dat.checkradius =true;
										dat.temproraryparticle = dat.intervals[interval].usingparticle;
										dat.particles.Particles[dat.temproraryparticle]->Position(dat.faces[usingface].constrictioncentre);		// put particle at the centre of constriction
										checkdistance (dat, dat.temproraryparticle);
										checkoverlap(dat, dat.temproraryparticle);
										if ((dat.checkdistance) and (!dat.checkoverlap))
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
						//if ((dat.checkdistance)and(dat.checkradius)and(!dat.checkoverlap))
						//	{
						//		dat.checkoverlap =true;
						//		for (int i=0; i<3; i++)
						//			{
						//				if (dat.overlappingpoint == dat.faces[m0].points[i])
						//					{
						//						cout << "hi \n";
						//						dat.checkoverlap =false;
						//					}
						//			}
						//		if (dat.checkoverlap)
						//		{
						//			dat.temproraryparticle= dat.overlappingpoint;
						//			dat.particlereuse = true;
						//		}
						//	}
					}
			}
	}
inline void useparticle(data & dat, int p0)
	{
		dat.particlesuse[p0]=true;
		dat.numberunusedparticles -=1;
		int interval = dat.particles.Particles[p0]->Tag;
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
