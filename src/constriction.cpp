#include <mechsys/dem/domain.h>
#include <iostream>
#include <cmath>
// MechSys
#include <mechsys/mesh/delaunay.h>
using namespace std;
using namespace DEM;
using namespace Mesh;
struct face
{
	int points[3];
	double constriction[3];
	Vec3_t constrictioncentres[3];
	int overlappoints[3];
	bool divided;
};
struct UserData
{
	DEM::Domain specimen;
	vector <face> faces;
	int count;
	int count1;
	int overlappingpoints[3];
	vector <int> specimen2mesh;
	vector <int> mesh2specimen;
};

int checkoverlap(void *UD, double constriction, Vec3_t constrictioncentre, int numberprocessors=6);
void * check (void *UD);
void matchid(void *UD, Mesh::Delaunay&mesh);
void planeconstriction (void *UD, int p[3], double&constriction, Vec3_t&constrictioncentre, int p0, Vec3_t x, Vec3_t Y);
void reduceconstriction (void * UD, int p[3], int overlapID, Vec3_t X, Vec3_t Y);
void sphereconstriction (int points[3], void * UD);
void textout(int p[3], void * UD);
int main(int argc, char **argv) try
{
	string filename;
	if (argc<2)
		{
			filename = "CSD";
		}
	else
		{
			filename =argv[1];

		}
	UserData ud;
	ud.specimen.Load(filename.c_str());
	Array <double> X(ud.specimen.Particles.Size());
	Array <double> Y(ud.specimen.Particles.Size());
	Array <double> Z(ud.specimen.Particles.Size());
	for (size_t i=0; i <ud.specimen.Particles.size();i++)
	{
		X[i]=ud.specimen.Particles[i]->x(0);
		Y[i]=ud.specimen.Particles[i]->x(1);
		Z[i]=ud.specimen.Particles[i]->x(2);
	}
	Mesh::Delaunay mesh(3);
	mesh.AddPoints(X,Y,Z);
	mesh.Generate();
	mesh.WriteVTU("CSD");
	size_t nc = mesh.Cells.Size();
	int p[3];
	ud.count=0;
	ud.count1=0;
//	cout << mesh.Cells[100]->V[1]->ID << " ";
//	for (int i=0; i<3; i++)
//	{
//		cout << mesh.Cells[100]->V[1]->C[i]<< " ";
//	}
//	cout << "\n";
//	cout << ud.specimen.Particles[mesh.Cells[100]->V[1]->ID]->x;
	matchid(&ud,mesh);
	for (size_t i=0;i<nc;i++)
	{
		for (int j=0;j<3;j++)
		{
			p[j]=ud.mesh2specimen[mesh.Cells[i]->V[j]->ID];
			cout << mesh.Cells[i]->V[j]->C[0] << "\n";
			cout << ud.specimen.Particles[ud.mesh2specimen[mesh.Cells[i]->V[j]->ID]]->x(0);
		}
		textout (p,&ud);
		sphereconstriction(p,&ud);
	}
	cout << "constriction " <<ud.faces.size()<<"\n";
	cout << "overlap " << ud.count <<"\n";
	cout << "no constriction " << ud.count1 << "\n";
}
MECHSYS_CATCH

inline int checkoverlap(void * UD, double constriction, Vec3_t constrictioncentre, int numberprocessors)
{
	UserData & ud = (*static_cast<UserData*>(UD));
	bool overlapcheck=false;
	int overlap;
	double overlapdistance = 0.;
	double distance = 0.;
	for (int i =0; i<ud.specimen.Particles.Size(); i++)
		{
			distance = norm(ud.specimen.Particles[i]->x-constrictioncentre) - constriction - ud.specimen.Particles[i]->Props.R;
			if (distance < overlapdistance)
			{
				overlapcheck=true;
				overlap =i;
				overlapdistance =distance;
			}
		}
	if (!overlapcheck)
	{
		overlap =-1;
	}
	return(overlap);
//	pthread_t thds[numberprocessors];
//	int rc;
//	for (long i=0;i<numberprocessors;i++)
//	{
//		rc = pthread_create(&thds[i], NULL, check, (void *)i);
//	}
}
inline void matchid (void *UD, Mesh::Delaunay&mesh)
{
	UserData & ud = (*static_cast<UserData*>(UD));
	if (!(mesh.Verts.Size()==ud.specimen.Particles.Size())) throw new Fatal("Number of vertexes and particles are not equal");
	bool check;
	for (size_t i=0; i<mesh.Verts.Size(); i++)
	{
		ud.mesh2specimen.push_back(-1);
		ud.specimen2mesh.push_back(-1);
	}
	for (size_t i=0; i<mesh.Verts.Size(); i++)
	{
		check =false;
		for (size_t j=0; j<ud.specimen.Particles.Size();j++)
		{
			if (mesh.Verts[i]->C[0]==ud.specimen.Particles[j]->x(0))
			{
				if ((mesh.Verts[i]->C[1]==ud.specimen.Particles[j]->x(1))and(mesh.Verts[i]->C[2]==ud.specimen.Particles[j]->x(2)))
				{
					check=true;
					ud.mesh2specimen[i]=j;
					ud.specimen2mesh[j]=i;
				}
			}
		}
//		if (!check)
//		{
//			cout << "Vertex " << i <<" not match \n";
//		}
	}
}
inline void planeconstriction(void *UD, int p[3], double&constriction, Vec3_t&constrictioncentre, int p0, Vec3_t X, Vec3_t Y)
{
	UserData & ud = (*static_cast<UserData*>(UD));
	double x[3];
	double y[3];
	double r[3];
	for (int i=0; i<3;i++)
	{
		x[i] = dot(ud.specimen.Particles[p[i]]->x-ud.specimen.Particles[p0]->x,X);
		y[i] = dot(ud.specimen.Particles[p[i]]->x-ud.specimen.Particles[p0]->x,Y);
		r[i] = ud.specimen.Particles[p[i]]->Props.R;
	}
	double z3 = pow(pow(norm(ud.specimen.Particles[p[2]]->x-ud.specimen.Particles[p0]->x),2.)-x[2]*x[2]-y[2]*y[2],0.5);
	double a1 = (r[1]*r[1]-r[0]*r[0]+x[0]*x[0]+y[0]*y[0]-x[1]*x[1]-y[1]*y[1])/2;
	double a2 = (r[2]*r[2]-r[0]*r[0]+x[0]*x[0]+y[0]*y[0]-x[2]*x[2]-y[2]*y[2]-z3*z3)/2;
	double b1 = ((r[2]-r[0])*(x[0]-x[1])-(r[1]-r[0])*(x[0]-x[2]))/((y[0]-y[2])*(x[0]-x[1])-(y[0]-y[1])*(x[0]-x[2]));
	double b2 = (a2*(x[0]-x[1])-a1*(x[0]-x[2]))/((y[0]-y[2])*(x[0]-x[1])-(y[0]-y[1])*(x[0]-x[2]));
	double c1 = (r[1]-r[0]-(y[0]-y[1])*b1)/(x[0]-x[1]);
	double c2 = (a1-(y[0]-y[1])*b2)/(x[0]-x[1]);
	double a = c1*c1+b1*b1-1;
	double b = c1*(c2-x[0])+b1*(b2-y[0])-r[0];
	double c = pow(c2-x[0],2.)+pow(b2-y[0],2.)-r[0]*r[0];
	if (b*b-a*c >0)
	{
		constriction = (-b+pow(b*b-a*c,0.5))/a;
		if (constriction>0)
		{
			constrictioncentre = X*(c1*constriction+c2)+Y*((b1*constriction+b2));//root of coordinate
		}
		else
		{
			constriction = (-b-pow(b*b-a*c,0.5))/a;
			if (constriction >0)
			{
				constrictioncentre = X*(c1*constriction+c2)+Y*((b1*constriction+b2));//root of coordinate
				cout << "Swap plane constriction \n";
			}
			else
			{
				constriction =-1;
				cout << "Plane constriction <0 \n";
			}
		}
	}
	else
	{
		cout << "no plane constriction \n";
		constriction =-1;
	}
}
inline void reduceconstriction( void * UD, int p[3], int overlapID, Vec3_t X, Vec3_t Y)
{
	UserData & ud = (*static_cast<UserData*>(UD));
	double constriction;
	Vec3_t constrictioncentre;
	face tempface;
	tempface.divided =true;
	int p1[3];
	for (int i=0; i<3; i++)
	{
		p1[0]=p[i];
		p1[1]=p[(i+1)%3];
		int checkID =overlapID;
		int count =0;
		while (checkID>-1)
		{
			p1[2]=checkID;
			planeconstriction(&ud, p1,constriction,constrictioncentre, p[0],X,Y);
			if (constriction>0)
			{
				checkID=checkoverlap(&ud, constriction,constrictioncentre);
			}
			else
			{
				break;
			}
			count+=1;
			if (count >10)
			{
				cout << " endless circle \n";
				break;
			}
		}
		if (constriction >0)
		{
			tempface.constriction[i]=constriction;
			tempface.constrictioncentres[i]=constrictioncentre;
			tempface.overlappoints[i]=p1[2];
		}
		else
		{
			tempface.constriction[i]=-1;
		}
	}
}
inline void sphereconstriction (int points[3], void * UD)
{
	UserData & ud = (*static_cast<UserData*>(UD));
	Vec3_t constrictioncentre;
	double r[3];
	double constriction;
	int p[3];
	for (int i =0;i<3;i++)
	{
		p[i]=points[i];
	}

	int swap;
	for (int i=0;i<2;i++)
	{
		for (int j=i+1;j<3;j++)
		{
			if (ud.specimen.Particles[p[i]]->Props.R<ud.specimen.Particles[p[j]]->Props.R)
			{
				swap = p[i];
				p[i]=p[j];
				p[j]=swap;
			}
		}
	}
	for (int i=0;i<3;i++)
	{
		r[i]=ud.specimen.Particles[p[i]]->Props.R;
	}
	Vec3_t X = ud.specimen.Particles[p[1]]->x - ud.specimen.Particles[p[0]]->x;
	double x2 = norm(X);
	X/=x2;
	double x3 = dot(ud.specimen.Particles[p[2]]->x - ud.specimen.Particles[p[0]]->x, X);
	Vec3_t Y = ud.specimen.Particles[p[2]]->x - ud.specimen.Particles[p[0]]->x-x3*X;
	double y3 = norm(Y);
	Y/=y3;
	double a1 = (r[0]-r[1])/x2;
	double b1 = (x2*x2+r[0]*r[0]+r[1]*r[1])/2/x2;
	double a2 = (r[0]-r[2]-a1*x3)/y3;
	double b2 = (x3*x3+y3*y3+r[0]*r[0]-r[2]*r[2]-2*b1*x3)/2/y3;
	double a = a1*a1+b1*b1-1;
	double b = a1*b1+a2*b2-r[0];
	double c = b1*b1+b2*b2-r[0]*r[0];
	if (b*b-a*c>=0)
	{
		constriction = (-b-pow(b*b-a*c,0.5))/a;
		if (constriction >0)
		{
			constrictioncentre=(a1*constriction+b1)*X+(a2*constriction+b2)*Y;
			textout (p,&ud);
			cout << "radius "<< constriction << " centre " << constrictioncentre << "\n";
		}
		else
		{
			constriction= (-b+pow(b*b-a*c,0.5))/a;
			if (constriction >0)
			{
				constrictioncentre=(a1*constriction+b1)*X+(a2*constriction+b2)*Y;
				cout << "swap constriction \n";
			}
			else
			{
				cout << "constriction <=0 \n";
				//textout(p,&ud);
				constriction =-1;
			}
		}
	}
	else
	{
		constriction =-1;
		cout << "no constriction \n";
		//textout(p, &ud);
	}
//	if (constriction >0)
//	{
//		int checkid;
//		checkid = checkoverlap(&ud,constriction, constrictioncentre);
//		if (checkid>-1)
//		{
//			reduceconstriction(&ud, p, checkid, X, Y);
//		}
//		else
//		{
//			face tempface;
//			tempface.constriction[0]= constriction;
//			tempface.constrictioncentres[0] = constrictioncentre;
//			for (int i=0; i<3; i++)
//			{
//				tempface.points[i]=p[i];
//			}
//			tempface.divided =false;
//			ud.faces.push_back(tempface);
//		}
//
//	}
}
inline void textout(int p[3], void * UD)
{
	UserData & ud = (*static_cast<UserData*>(UD));
	for( int i =0; i<3;i++)
	{
		cout << "radius: " << ud.specimen.Particles[p[i]]->Props.R << " centre: " << ud.specimen.Particles[p[i]]->x <<"\n";
	}
	int hi;
	cin >> hi;
}
//inline void * check (void *UD)
//{
//
//}

