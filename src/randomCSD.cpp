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
	double x1;
	double x2;
	double y2;
	Vec3_t normal;
	Vec3_t centre;
	Vec3_t constrictioncentre;
	Vec3_t systemcoordinate[3];


};
struct data
{
	bool choice;
	bool check;
	bool checkin;
	bool facereuse;
	bool particlereuse;
	bool checkoverlap;
	bool checkdistance;
	bool checkradius;
	bool checkboundary;
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
	int mininterval;
	int maxinterval;
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
	face facetemprorary;
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
void observeface(data&dat, int m0);
void putparticle(data & dat, int p0, int m0);
void screenreport(data&dat);
void savedomain(data&dat);
void textout (data& dat);
void tryputparticle(data&dat, int p0,int m0);
void useparticle(data & dat, int p0);
double interpolation (double x1,double x2, double y1,double y2,double x3);
int randomselect(data&dat, int max);
double constriction(data&dat, int m0);
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
	int usingparticle;
	cout << "Reading data \n";
	ifstream datain;
	//cout << "Welcome to CSD program \n";
	//cout << "Enter data file name: ";
	//cin >> dat.datafile;
	dat.datafile = "sand.csd";
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
	cout << dat.maxfraction <<"\n";
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
	cout << "Calculating intervals \n";
	dat.choice = false;
	while (!dat.choice)
		{
			calculateintervals(dat);																														// generate list of intervals with number of particles
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
	cout << "preparing list \n";
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
	// generate basic tetrahedron
	cout << "generating basic tetrahedron \n";
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
	dat.numberfaces =0;
	dat.numberopenfaces =0;
	createface(dat, dat.usingparticle, dat.usingparticle-1, dat.usingparticle-2,dat.usingparticle-3);								//create face from 3 particles
	dat.faces[0].faceuse=-1;
	putparticle(dat,dat.usingparticle-3,0);																							//put first particle on first face, this face will put on twice
	createtetrahedron(dat,0,dat.usingparticle-3);
	dat.faces[0].faceuse=1;																											// mark used side
	dat.usingface=0;

	for (int i=0; i<4; i++)																							//mark used particles
		{
			useparticle(dat, dat.usingparticle);
		}
	// first loop for putting particles
	cout << "first loop running \n";
	int usinginterval = dat.usinginterval;
	while ((dat.numberunusedparticles> 0)and(dat.numberopenfaces > 0))
		{
			while (dat.particlesuse[dat.usingparticle])
				{
					dat.usingparticle -=1;
					if (dat.usingparticle==0)
						{
							cout << "using particle =0 \n";
							break;
						}
				}
			while (dat.intervals[dat.usinginterval].ability ==false)
				{
					dat.usinginterval -=1;
					if (dat.usinginterval ==0)
						{
							cout << "using interval =0 \n";
							break;
						}
				}
			while (dat.faces[dat.usingface].faceuse == 0)
				{
					dat.usingface +=1;
					if (dat.usingface == dat.numberfaces-1)
						{
							cout << "using face = number faces -1 \n";
							screenreport(dat);
							break;
						}
				}
			//cout << dat.numberunusedparticles <<"\n";
			if (dat.numberunusedparticles ==17144)
				{
					dat.check=true;
				}
			usingparticle =randomselect(dat, dat.usinginterval);														// select particle randomly
			tryputparticle(dat,usingparticle, dat.usingface);
			closeface(dat, dat.usingface);
			// build new tetrahedron with new faces
			if (dat.checkradius and dat.checkdistance)
				{
					createtetrahedron(dat, dat.usingface, dat.temproraryparticle);										// create new tetrahedron
					//screenreport(dat);
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

		}
	for (int i=0; i< dat.particles.Particles.Size(); i++)
		{
			if (dat.particlesuse[i]==false)
			{
				dat.particles.Particles[i]->Tag -=100;
			}

		}
	screenreport(dat);
	// export results
	savedomain(dat);
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
inline int randomselect(data&dat, int max)
	{
		int interval;
		int rs;
		do
			{
				interval =rand()%(dat.maxinterval-dat.mininterval)+dat.mininterval;
			}
		while (dat.intervals[interval].ability ==false);
		return (dat.intervals[interval].usingparticle);
	}
inline double constriction(data&dat, int m0)
	{
		dat.faces[m0].x1 = norm(dat.particles.Particles[dat.faces[m0].points[1]]->x-dat.particles.Particles[dat.faces[m0].points[0]]->x);	// get coordinates of particles in local coordinate system
		dat.faces[m0].x2 = dot(dat.particles.Particles[dat.faces[m0].points[2]]->x - dat.particles.Particles[dat.faces[m0].points[0]]->x, dat.faces[m0].systemcoordinate[0]);
		dat.faces[m0].y2 = dot(dat.particles.Particles[dat.faces[m0].points[2]]->x - dat.particles.Particles[dat.faces[m0].points[0]]->x, dat.faces[m0].systemcoordinate[1]);
		double a1 = (dat.particles.Particles[dat.faces[m0].points[0]]->Props.R-dat.particles.Particles[dat.faces[m0].points[1]]->Props.R)/dat.faces[m0].x1;
		double a2 = (dat.particles.Particles[dat.faces[m0].points[0]]->Props.R -dat.particles.Particles[dat.faces[m0].points[2]]->Props.R-a1*dat.faces[m0].x2)/dat.faces[m0].y2;
		double b1 = (pow(dat.faces[m0].x1,2.0)+pow(dat.particles.Particles[dat.faces[m0].points[0]]->Props.R,2.0)-pow(dat.particles.Particles[dat.faces[m0].points[1]]->Props.R,2.0))/2.0/dat.faces[m0].x1;
		double b2 = (pow(dat.faces[m0].x2,2.0)+pow(dat.faces[m0].y2,2.0)+pow(dat.particles.Particles[dat.faces[m0].points[0]]->Props.R,2.0)-pow(dat.particles.Particles[dat.faces[m0].points[2]]->Props.R,2.0)-2*b1*dat.faces[m0].x2)/2/dat.faces[m0].y2;
		double a = pow(a1,2.0)+pow(a2,2.0)-1;
		double b = a1*b1+a2*b2-dat.particles.Particles[dat.faces[m0].points[0]]->Props.R;
		double c = pow(b1, 2.0) +pow(b2, 2.0) - pow(dat.particles.Particles[dat.faces[m0].points[0]]->Props.R, 2.0);
		if ((b*b-a*c)<0)
			{
				//cout << "Constriction problems \n";
				//observeface(dat,m0);
				dat.position = (dat.particles.Particles[dat.faces[m0].points[0]]->x+dat.particles.Particles[dat.faces[m0].points[1]]->x+dat.particles.Particles[dat.faces[m0].points[3]]->x)/3;
				return (100.0);
			}
		double constriction = (-b-pow(pow(b, 2.0)-a*c,0.5))/a;
		if (constriction >0)
			{
				dat.position = (a1*constriction+b1)*dat.faces[m0].systemcoordinate[0]+(a2*constriction+b2)*dat.faces[m0].systemcoordinate[1]+dat.particles.Particles[dat.faces[m0].points[0]]->x;
				return (constriction);
			}
		else
			{
				constriction = (-b+pow(pow(b, 2.0)-a*c,0.5))/a;
				if (constriction >0)
					{
						dat.position = (a1*constriction+b1)*dat.faces[m0].systemcoordinate[0]+(a2*constriction+b2)*dat.faces[m0].systemcoordinate[1]+dat.particles.Particles[dat.faces[m0].points[0]]->x;
						return (constriction);
					}
				else
					{
						//cout << " both constriction size are smaller than 0 \n";
						//observeface(dat,m0);
						dat.position = (dat.particles.Particles[dat.faces[m0].points[0]]->x+dat.particles.Particles[dat.faces[m0].points[1]]->x+dat.particles.Particles[dat.faces[m0].points[3]]->x)/3;
						return (100.0);
					}
			}
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
		dat.checkboundary =false;
		if (dat.specimentype.compare("Cube")==0)
			{
				for (int i = 0; i < 3; i++)
					{
						if (abs(dat.particles.Particles[p0]->x(i))+dat.particles.Particles[p0]->Props.R > dat.specimen[i])
							{
								dat.checkdistance= false;
							}
						else
							{
								dat.checkboundary =true;
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
		Array<Vec3_t> V(3);
		for (int i=0;i<3;i++)
			{
				dat.facetemprorary.points[i]=sort[i];
				V[i]=dat.particles.Particles[sort[i]]->x;
			}
		if (dat.particlereuse)																																	// check whether face is reused or not
			{
				for (int i=0; i< dat.numberfaces; i++)
					{
						if (dat.facetemprorary.points[0]== dat.faces[i].points[0])
							{
								if ((dat.facetemprorary.points[1]== dat.faces[i].points[1])and(dat.facetemprorary.points[2]== dat.faces[i].points[2]))
									{
										dat.facereuse=true;																										// *********does it need to calcualte facetemprorary.faceuse to know which side particle is put ?
										dat.facereusenumber = i;
									}
							}
					}
			}
		if (!dat.facereuse)
			{
				dat.faces.push_back(dat.facetemprorary);																										// add user-defined face
				dat.numberfaces +=1;
				dat.numberopenfaces +=1;
				dat.faces[dat.numberfaces-1].point = p3;
				DEM::Face facedemtemprorary(V);																												// add face in library
				facedemtemprorary.Normal(dat.faces[dat.numberfaces-1].normal);
				establishlocalsystem(dat, dat.numberfaces-1);
				for (int i=0; i<3; i++)
					{
						dat.faces[dat.numberfaces-1].systemcoordinate[i]=dat.localsystem[i];
					}
				dat.faces[dat.numberfaces-1].constrictionsize = constriction(dat,dat.numberfaces-1);																			// calculate constriction size
				dat.faces[dat.numberfaces-1].constrictioncentre = dat.position;
				dat.faces[dat.numberfaces-1].normal /= norm(dat.faces[dat.numberfaces-1].normal);
				//observeface(dat, dat.numberfaces-1);
				double distance = dot(dat.particles.Particles[dat.faces[dat.numberfaces-1].point]->x - dat.particles.Particles[dat.faces[dat.numberfaces-1].points[0]]->x, dat.faces[dat.numberfaces-1].normal);
				if (distance>0)
					{
						dat.faces[dat.numberfaces-1].faceuse = 1;
					}
				else
					{
						dat.faces[dat.numberfaces-1].faceuse =-1;
					}
			}
		else
			{
				dat.facetemprorary.point = p3;
				DEM::Face facedemtemprorary(V);																												// add face in library
				facedemtemprorary.Normal(dat.facetemprorary.normal);
				double distance = dot(dat.particles.Particles[dat.facetemprorary.point]->x - dat.particles.Particles[dat.facetemprorary.points[0]]->x, dat.facetemprorary.normal);
				if (distance>0)
					{
						dat.facetemprorary.faceuse = 1;
					}
				else
					{
						dat.facetemprorary.faceuse =-1;
					}
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
				if (dat.facetemprorary.faceuse*dat.faces[dat.facereusenumber].faceuse <0)
					{
						closeface(dat,dat.facereusenumber);
					}
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
				if (dat.facetemprorary.faceuse*dat.faces[dat.facereusenumber].faceuse <0)
					{
						closeface(dat,dat.facereusenumber);
					}
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
				if (dat.facetemprorary.faceuse*dat.faces[dat.facereusenumber].faceuse <0)
					{
						closeface(dat,dat.facereusenumber);
					}
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
		dat.mininterval = 0;
		dat.maxinterval = dat.numberintervals -1;
	}
inline void moveon(data & dat, int p0, int m0)
	{
		Vec3_t distance = dat.particles.Particles[p0]->x - dat.particles.Particles[dat.overlappingpoint]->x;
		Vec3_t projection = dot(distance, dat.faces[m0].normal);
		double move = sqrt(pow(dat.particles.Particles[p0]->Props.R + dat.particles.Particles[dat.overlappingpoint]->Props.R,2.0) - pow(norm(distance),2.0)+pow(norm(projection),2.0))-norm(projection);
		dat.finalposition = dat.particles.Particles[p0]->x - dat.faces[m0].faceuse*move*dat.faces[m0].normal;
		dat.particles.Particles[p0]->Position(dat.finalposition);
	}
inline void observeface(data&dat, int m0)
	{
		int p1 = dat.faces[m0].points[0];
		int p2 = dat.faces[m0].points[1];
		int p3 = dat.faces[m0].points[2];
		int p4 = dat.faces[m0].point;
		Vec3_t e1 = dat.particles.Particles[dat.faces[m0].points[0]]->x;
		Vec3_t e2 = dat.particles.Particles[dat.faces[m0].points[1]]->x;
		Vec3_t e3 = dat.particles.Particles[dat.faces[m0].points[2]]->x;
		Vec3_t e4 = dat.particles.Particles[dat.faces[m0].point]->x;
		double r1 = dat.particles.Particles[dat.faces[m0].points[0]]->Props.R;
		double r2 = dat.particles.Particles[dat.faces[m0].points[1]]->Props.R;
		double r3 = dat.particles.Particles[dat.faces[m0].points[2]]->Props.R;
		double r4 = dat.particles.Particles[dat.faces[m0].point]->Props.R;
		double x1 = dat.faces[m0].x1;
		double x2 = dat.faces[m0].x2;
		double y2 = dat.faces[m0].y2;
		double constriction = dat.faces[m0].constrictionsize;
		Vec3_t normal = dat.faces[m0].normal;

		dat.check= true;
	}
inline void screenreport(data&dat)
	{
		cout << "Number of particles: " << dat.numberparticles<< "\n ";
		cout << "  		Number of unused particles: " << dat.numberunusedparticles << "\n";
		cout << "  		Current using particle: "<< dat.usingparticle <<"\n ";
		cout << "  	Number of faces: "<<dat.numberfaces<< "\n";
		cout << "  		Current using face: "<< dat.usingface<< "\n";
		cout << "  		Number of open faces: "<<dat.numberopenfaces<<"\n";
		dat.check =true;
	}
inline void putparticle(data & dat, int p0, int m0)
	{
		dat.checkradius =false;
		double r[4];
		for (int i =0;i<3;i++)
			{
				r[i] = dat.particles.Particles[dat.faces[m0].points[i]]->Props.R;
			}
		r[3]=dat.particles.Particles[p0]->Props.R;
		double a1 = (r[0] -r[1])/dat.faces[m0].x1;
		double b1 = (pow(dat.faces[m0].x1,2.)+pow(r[0],2.)-pow(r[1],2.))/2./dat.faces[m0].x1;
		double a2 = (r[0]-r[2]-a1*dat.faces[m0].x2)/dat.faces[m0].y2;
		double b2 = (pow(dat.faces[m0].x2,2.)+pow(dat.faces[m0].y2,2.)+pow(r[0],2.)-pow(r[2],2.)-2.*b1*dat.faces[m0].x2)/2./dat.faces[m0].y2;
		double x4 = a1*r[3]+b1;
		double y4 = a2*r[3]+b2;
		double z4 = pow(r[0]+r[3],2)-pow(a1*r[3]+b1,2)-pow(a2*r[3]+b2,2);
		if (z4 > 0)
			{
				z4=-dat.faces[m0].faceuse*pow(z4,0.5);
				dat.particles.Particles[p0]->Position(x4*dat.faces[m0].systemcoordinate[0]+y4*dat.faces[m0].systemcoordinate[1]+z4*dat.faces[m0].systemcoordinate[2]);
				dat.checkradius =true;
			}
		//else
		//	{
		//		cout << "putting problems \n";
		//	}
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
		if (dat.numberopenfaces <10)
			{
				dat.check=true;
			}
		dat.temproraryparticle = p0;
		int usingface = m0;
		double distancepoint =0;
		double distanceoverlap =0;																									// distance from overlapping point
		putparticle(dat,dat.temproraryparticle,usingface);																			// put particle first time
		checkoverlap(dat, dat.temproraryparticle);
		checkdistance(dat, dat.temproraryparticle);
		if ((!dat.checkradius)or(!dat.checkdistance)or(!dat.checkoverlap))
			{
				int interval = dat.particles.Particles[dat.temproraryparticle]->Tag;
				//distancepoint = dot(dat.particles.Particles[dat.faces[usingface].point]->x - dat.particles.Particles[dat.faces[usingface].points[0]]->x, dat.faces[usingface].normal);	//distance from vertice of tetrahedron
				if (dat.particles.Particles[dat.temproraryparticle]->Props.R > dat.faces[usingface].constrictionsize)				// if particle size is bigger than constriction size
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
										if ( (!dat.checkdistance) or (!dat.checkoverlap))
											{
												distanceoverlap = dot(dat.particles.Particles[dat.overlappingpoint]->x - dat.particles.Particles[dat.faces[usingface].points[0]]->x, dat.faces[usingface].normal);
												while ((dat.faces[usingface].faceuse*distanceoverlap> -dat.approximation)and(!dat.checkoverlap)and(dat.checkdistance))							// if overlap with particle from other side of face then move particle
													{
														moveon(dat,dat.temproraryparticle,usingface);
														checkoverlap(dat, dat.temproraryparticle);
														checkdistance(dat,dat.temproraryparticle);
														if (!dat.checkoverlap) //how about not overlap anymore ?
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
										double x =0;
										double y=0;
										if (dat.faces[m0].x1-dat.particles.Particles[dat.faces[m0].points[0]]->Props.R - dat.particles.Particles[dat.faces[m0].points[1]]->Props.R -dat.particles.Particles[dat.temproraryparticle]->Props.R<-dat.approximation) // if p0 and p1 are close enough
											{
												x= (pow(dat.particles.Particles[dat.faces[m0].points[0]]->Props.R+dat.particles.Particles[dat.temproraryparticle]->Props.R,2.0)-pow(dat.particles.Particles[dat.faces[m0].points[1]]->Props.R+dat.particles.Particles[dat.temproraryparticle]->Props.R,2.0)+pow(dat.faces[m0].x1,2.0))/2.0/dat.faces[m0].x1;
												y= pow(pow(dat.particles.Particles[dat.faces[m0].points[0]]->Props.R+dat.particles.Particles[dat.temproraryparticle]->Props.R,2.0)-pow(x,2.0),0.5);
												dat.particles.Particles[dat.temproraryparticle]->x = x*dat.faces[m0].systemcoordinate[0]+y*dat.faces[m0].systemcoordinate[1]+dat.particles.Particles[dat.faces[m0].points[0]]->x;				// add particles on the face
											}
										else
											{
												dat.particles.Particles[dat.temproraryparticle]->x = dat.particles.Particles[dat.faces[m0].points[0]]->x+dat.faces[m0].x1/2*dat.faces[m0].systemcoordinate[0]-2*dat.approximation*dat.faces[m0].faceuse*dat.faces[m0].systemcoordinate[3];
											}
										checkdistance (dat, dat.temproraryparticle);
										checkoverlap(dat, dat.temproraryparticle);
										if ((!dat.checkdistance) or (!dat.checkoverlap))
											{
												distanceoverlap = dot(dat.particles.Particles[dat.overlappingpoint]->x - dat.particles.Particles[dat.faces[usingface].points[0]]->x, dat.faces[usingface].normal);
												while ((dat.faces[usingface].faceuse*distanceoverlap> -dat.approximation)and (!dat.checkoverlap)and(dat.checkdistance))							// if overlap with particle from other side of face then move particle
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
						if ((dat.checkdistance)and(dat.checkradius)and(!dat.checkoverlap))
							{
								dat.temproraryparticle= dat.overlappingpoint;
								dat.particlereuse = true;
							}
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
		if (p0== dat.usingparticle)
			{
				dat.usingparticle -=1;
			}
		while ((dat.intervals[dat.mininterval].ability ==false)and (dat.mininterval<dat.maxinterval))
			{
				dat.mininterval+=1;
			}
		while ((dat.intervals[dat.maxinterval].ability ==false)and (dat.mininterval<dat.maxinterval))
			{
				dat.maxinterval-=1;
			}
	}
