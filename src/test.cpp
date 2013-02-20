#include <mechsys/dem/domain.h>
// variable declaration
using namespace std;
using namespace DEM;
int main()
{
	DEM::Domain test;
	test.AddSphere(-1,Vec3_t (0.0,0.0,0.0), 1.0, 3.2);
	test.AddSphere(-2,Vec3_t (2.0,0.0,0.0),1.0,3.2);
	test.AddSphere(-2,Vec3_t (1.0,pow(3.0,0.5),0.0),1.0,3.2);
	test.AddSphere(-2,Vec3_t (1.0,pow(3.0,-0.5),pow(8/3.0,0.5)),1.0,3.2);
	test.WriteXDMF("test");

}
