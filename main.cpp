#include <iostream>
#include "consts.h"
#include "func.h"
#include "classes.h"
#include "vtk_cpp.h"

#define imax 200
#define jmax 200
#define nmax 5000

const double dx = 1.0/143.0;
const double dt = 0.5;

using namespace std;

int main()
{
 GS *gs = new GS [imax*jmax];

 initialize(imax,jmax,gs);

 for (int n = 0; n<nmax; n++)
 {
  cout<<"n = "<<n<<endl;
  advance(imax,jmax,dx,dt,gs);
 }

 output_vtk(imax,jmax,dx,gs);

 delete [] gs;
}
