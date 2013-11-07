#include <mechsys/dem/domain.h>
// variable declaration
using namespace std;
using namespace DEM;
int main(int argc, char **argv) try
{
	string filename;
	DEM::Domain particles;
	vector <Vec3_t> centers;
	vector <double> radii;
	Vec3_t center;
	double radius;
	double a[4];
	double distance;
	double divide;
	string domainin;
	string domainout;
	ifstream datain;
	ofstream dataout;
	if (argc<2)
		{
			filename = "slicesphere.drp";
		}
	else
		{
			filename = argv[1];
		}

	datain.open(filename.c_str());
	datain >> domainin; 				datain.ignore(200,'\n');
	datain >> domainout;				datain.ignore(200,'\n');
//	domainout.Printf    ("%s.%s", domainout.c_str(), "out");
	for (int i=0; i<4; i++)
		{
			datain >> a[i];
		}
										datain.ignore(200,'\n');
	datain.close();
	divide = pow(a[0]*a[0]+a[1]*a[1]+a[2]*a[2],0.5);
	for (int i=0;i<4;i++)
		{
			a[i]/=divide;
		}
	Vec3_t normvec (a[0],a[1],a[2]);
	particles.Load(domainin.c_str());
	dataout.open(domainout.c_str());
	for (size_t i=0; i<particles.Particles.Size();i++)
		{
			distance=a[3]+dot(normvec,particles.Particles[i]->x);
			if (abs(distance) <particles.Particles[i]->Props.R)
				{
					center = particles.Particles[i]->x -distance*normvec;
//					center.push_back(center);
//					radii.push_back(pow(pow(particles.Particles[i]->Props.R,2.)+distance*distance,0.5));
					radius =pow(pow(particles.Particles[i]->Props.R,2.)-pow(distance,2.),0.5);
					for (int j=0;j<3;j++)
						{
							dataout << center(j)<< " ";
						}
					dataout << radius << "\n";
				}
		}
	dataout.close();
}
MECHSYS_CATCH
