#include <cmath>
#include <iostream>
#include <limits>
#include <cstdlib>
#include <ctime>
#include "consts.h"
#include "classes.h"
using namespace std;

void initialize(int i_max, int j_max, GS *gs)
{
 int i, j, l, k, ic, jc, is, ie, js, je;
 double r1, r2, uu, vv;
 const int nislands = 35;
 const int rad1 = 5;

 for (i=0; i<i_max; i++)
 {
  for (j=0; j<j_max; j++)
  {
   l = j + i*j_max;
   gs[l].edit_u(0.42);
   gs[l].edit_v(0.28);
  }
 }

 for (k=0; k<nislands; k++)
 {
 
   // choose center of islands randomly

   r1=rand() / double(RAND_MAX);
   r2=rand() / double(RAND_MAX);
   ic = int(r1*((double) i_max));
   jc = int(r2*((double) j_max));

   is = ic - rad1;
   ie = ic + rad1;
   is = is>=0 ? is : 0;
   ie = ie<=i_max-1? ie : i_max-1;

   js = jc - rad1;
   je = jc + rad1;
   js = js>=0 ? js : 0;
   je = je<=j_max-1? je : j_max-1;

   r1=rand() / double(RAND_MAX);
   r2=rand() / double(RAND_MAX);
   uu = r1;
   vv = r2;

   uu = uu >= 0.0? uu : 0.0;
   uu = uu <= 1.0? uu : 1.0;
   vv = vv >= 0.0? vv : 0.0;
   vv = vv <= 1.0? vv : 1.0;

   for (i=is; i<=ie; i++)
   {
   for (j=js; j<=je; j++)
   {
    l = j + i*j_max; 

    gs[l].edit_u(uu);
    gs[l].edit_v(vv);
   }
   }
 }

}

//-------------------------------------------------

void advance(int i_max, int j_max, double dx, double dt, GS *gs)
{
 GS *dgsdt = new GS [i_max*j_max];
 
 int l;
 double d2udx2, d2udy2, d2vdx2, d2vdy2;
 double Del2u, Del2v, dudt, dvdt, uu, vv;

 for (int i=0; i<i_max; i++)
 {
  for (int j=0; j<j_max; j++)
  {
   l = j + i*j_max;

   // compute second derivative
   if(j>0 && j<j_max-1)
   {
    d2udx2 = (gs[l-1].get_u() - 2.0*gs[l].get_u() + gs[l+1].get_u())/(dx*dx);   
    d2vdx2 = (gs[l-1].get_v() - 2.0*gs[l].get_v() + gs[l+1].get_v())/(dx*dx);   
   }     
   else if(j==0)
   {
    d2udx2 = (gs[l+j_max-1].get_u() - 2.0*gs[l].get_u() + gs[l+1].get_u())/(dx*dx);
    d2vdx2 = (gs[l+j_max-1].get_v() - 2.0*gs[l].get_v() + gs[l+1].get_v())/(dx*dx);
   }
   else if(j==j_max-1)
   {
    d2udx2 = (gs[l-1].get_u() - 2.0*gs[l].get_u() + gs[l-j_max+1].get_u())/(dx*dx);   
    d2vdx2 = (gs[l-1].get_v() - 2.0*gs[l].get_v() + gs[l-j_max+1].get_v())/(dx*dx);   
   }

   if(i>0 && i<i_max-1)
   {
    d2udy2 = (gs[l-j_max].get_u() - 2.0*gs[l].get_u() + gs[l+j_max].get_u())/(dx*dx);
    d2vdy2 = (gs[l-j_max].get_v() - 2.0*gs[l].get_v() + gs[l+j_max].get_v())/(dx*dx);
   }
   else if(i==0)
   {
    d2udy2 = (gs[l+(i_max-1)*j_max].get_u() - 2.0*gs[l].get_u() + gs[l+j_max].get_u())/(dx*dx);
    d2vdy2 = (gs[l+(i_max-1)*j_max].get_v() - 2.0*gs[l].get_v() + gs[l+j_max].get_v())/(dx*dx);
   }
   else if(i==i_max-1)
   {
    d2udy2 = (gs[l-j_max].get_u() - 2.0*gs[l].get_u() + gs[l-(i_max-1)*j_max].get_u())/(dx*dx);
    d2vdy2 = (gs[l-j_max].get_v() - 2.0*gs[l].get_v() + gs[l-(i_max-1)*j_max].get_v())/(dx*dx);
   }

   Del2u = d2udx2 + d2udy2;    
   Del2v = d2vdx2 + d2vdy2;    
   uu = gs[l].get_u(); 
   vv = gs[l].get_v(); 

   dudt = Du*Del2u - uu*vv*vv + F*(1.0-uu);
   dvdt = Dv*Del2v + uu*vv*vv - (F+k)*vv;
 
   dgsdt[l].edit_u(dudt); 
   dgsdt[l].edit_v(dvdt); 

  }
 } 

 
 for (int i=0; i<i_max; i++)
 {
  for (int j=0; j<j_max; j++)
  {
   l = j + i*j_max;

   uu = gs[l].get_u() + dgsdt[l].get_u()*dt;
   vv = gs[l].get_v() + dgsdt[l].get_v()*dt;

   gs[l].edit_u(uu);
   gs[l].edit_v(vv);
  }
 } 


 delete [] dgsdt; 
}

//-------------------------------------------------

