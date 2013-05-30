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
	int points[3];
	int point;
	int faceuse;
	double constrictionsize;
	Vec3_t normal;
	Vec3_t constriction;
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
void closeface(data&dat, int m0);
void convertcoordinates(data & dat, Vec3_t & position, bool way);
void createface(data & dat, int p0, int p1, int p2,int faceuse);
void createtetrahedron(data & dat, int m0, int p0);
void establishlocalsystem(data&dat, int m0);
void gsdgenerate(data & dat);
void listprepare(data & dat);
void moveon(data & dat, int p0);
void putparticle(data & dat, int p0, int m0);
void textout (data& dat);
void useparticle(data & dat, int p0);
double interpolation (double x1,double x2, double y1,double y2,double x3);
// main procedure
int main(int argc, char** argv) try
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
	int interval=0;
	int usingparticle =0;
	int usingface = 0;
	double distance1=0;
	double distance2=0;
	ifstream datain;
	//cout << "Welcome to CSD program \n";
	//cout << "Enter data file name: ";
	//cin >> dat.datafile;

	cout << "Reading data \n";
	cout << filename << "\n";
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
	datain >> dat.GSD;				datain.ignore(200,'\n');
	cout << dat. maxfraction << "\n";
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
	while (!dat.choice)
		{
			calculateintervals(dat);																														// generate list of intervals with number of particles
			cout << dat.numberparticles <<"\n";
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
	cout << "Create basic face \n";
	dat.numberunusedparticles = dat.numberparticles;
	// generate basic tetrahedron
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
	dat.numberunusedparticles -=4;
	usingparticle = dat.usingparticle;
	// first loop for putting particles
	cout << "Start first loop \n";
	int usinginterval = dat.usinginterval;
		while ((dat.usingparticle > -1)and(dat.numberopenfaces > 0))
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
				usingparticle = dat.usingparticle;
				usingface = dat.usingface;
				putparticle(dat,usingparticle,dat.usingface);
				checkoverlap(dat, usingparticle);
				checkdistance(dat, usingparticle);
				if (!dat.checkradius)
					{
						closeface(dat, dat.usingface);
					}
				else if (!dat.checkdistance)
					{
						if (dat.fillingtype.compare("Layer")==0)
							{
								closeface(dat, dat.usingface);												// ---------------------------should calculate pore space with boundary
							}
						else if (dat.fillingtype.compare("Random")==0)
							{
								closeface(dat, dat.usingface);
							}
					}
				else if (!dat.checkoverlap)
					{
						distance1 = dot(dat.particles.Particles[dat.faces[dat.usingface].point]->x - dat.particles.Particles[dat.faces[dat.usingface].points[0]]->x, dat.faces[dat.usingface].normal);
						while (!dat.checkoverlap)
							{
								distance2 = dot(dat.particles.Particles[dat.overlappingpoint]->x - dat.particles.Particles[dat.faces[dat.usingface].points[0]]->x, dat.faces[dat.usingface].normal);
								if (distance1*distance2<0)
									{
										usingparticle = dat.overlappingpoint;
										dat.particlereuse = true;
										dat.checkoverlap=true;
									}
								else
									{
										moveon(dat, usingparticle);
										checkoverlap(dat, usingparticle);
										checkdistance(dat, usingparticle);								// to be sure that moving particle is not out of boundary
										if (!dat.checkdistance)
											{
												closeface(dat, dat.usingface);
												dat.checkoverlap =true;									// now checkdistance is false
											}
									}
							}
					}
				// build new tetrahedron with new faces
				if (dat.checkradius and dat.checkoverlap and dat.checkdistance)
					{
						closeface(dat, dat.usingface);
						createtetrahedron(dat, dat.usingface, usingparticle);										// create new tetrahedron
						if (dat.particlereuse)
							{
								dat.particlereuse = false;
							}
						else
							{
								useparticle(dat, usingparticle);													// mark particle used
								dat.usingparticle -=1;
								dat.numberunusedparticles -=1;
							}
						//dat.intervals[usinginterval].usingparticle -=1;
						//if (dat.intervals[usinginterval].usingparticle < dat.intervals[usinginterval].firstparticle)
						//	{
						//		dat.intervals[usinginterval].ability=false;
						//		usinginterval -=1;
						//	}
					}
				dat.usingface+=1;
			}
		// second loop for putting particles
		//Mesh::Unstructured mesh(dat.particles.Particles.Size());
		//Array<double> X(dat.particles.Particles.Size()), Y(dat.particles.Particles.Size()), Z(dat.particles.Particles.Size());
		//for (size_t i=0;i<dat.particles.Particles.Size();i++)
		//	{
		//		X[i]=dat.particles.Particles[i]->x(0);
		//		Y[i]=dat.particles.Particles[i]->x(1);
		//		Z[i]=dat.particles.Particles[i]->x(2);
		//	}
		//mesh.Delaunay(X, Y, Z, -1);
		//cout << "saving mesh \n";
		//mesh.WriteVTU ("delaunay3d");

		// export results
		cout << "delete unused particles \n";
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
		cout << "save files \n";
		dat.particles.WriteXDMF(dat.outputdomain.c_str());// export to draw visual results
		dat.particles.Save(dat.outputdomain.c_str());
		textout (dat);
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
		else if (dat.specimentype.compare("Cylinder")==0)
			{
				if (abs(dat.particles.Particles[p0]->x(2))+dat.particles.Particles[p0]->Props.R > dat.specimen[1])
					{
						dat.checkdistance=false;
					}
				else if (abs(pow(pow(dat.particles.Particles[p0]->x(0),2.)+pow(dat.particles.Particles[p0]->x(1),2.),0.5))+dat.particles.Particles[p0]->Props.R>dat.specimen[0])
					{
						dat.checkdistance=false;
					}
			}
	}
inline void checkoverlap(data&dat, int p0)
	{
		dat.checkoverlap =true;
		dat.overlappingpoint =-1;
		double overlapdistance = -0.00001;																												// default overlappint distance, to avoid calculation error
		for ( int i = 0;  i < dat.numberparticles; i++)
			{
				if ((dat.particlesuse[i])and(!(i==p0)))
					{
						double distance = norm(dat.particles.Particles[p0]->x-dat.particles.Particles[i]->x) -dat.particles.Particles[p0]->Props.R -dat.particles.Particles[i]->Props.R;
						if (distance < - dat.approximation)
							{
								if (overlapdistance > distance)
									{
										dat.overlappingpoint = i;
										overlapdistance = distance;
									}
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
				dat.intervals[i].ability =true;\
				cout << i << " " << dat.intervals[i].numberparticles << "\n";
			}
	}
inline void moveon(data & dat, int p0)
	{
		Vec3_t distance = dat.particles.Particles[p0]->x - dat.particles.Particles[dat.overlappingpoint]->x;
		Vec3_t projection = dot(distance, dat.faces[dat.usingface].normal);
		double move = sqrt(pow(dat.particles.Particles[p0]->Props.R + dat.particles.Particles[dat.overlappingpoint]->Props.R,2) - pow(norm(distance),2)+pow(norm(projection),2))-norm(projection);
		dat.finalposition = dat.particles.Particles[p0]->x - dat.faces[dat.usingface].faceuse*move*dat.faces[dat.usingface].normal;
		dat.particles.Particles[p0]->Position(dat.finalposition);
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
inline void textout (data & dat)
	{
		ofstream dataout;																							// export to text file
		String fn;
		fn.Printf    ("%s.%s", dat.outputdomain.c_str(), "out");
		dataout.open(fn.c_str());
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
inline void useparticle(data & dat, int p0)
	{
		dat.particlesuse[p0]=true;
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
