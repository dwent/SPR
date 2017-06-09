/*************************************************************************
Improvements:
Inline Functions should place in head file. 2007.5.11.
BUG: PointsOnLine()
c3[0]=c2[0]-(c1[0]-c2[0])*dis32/dis12;=>c3[0]=c2[0]-(c1[0]-c2[0])*dis32/dis12;
Date 2008.7.11.

*************************************************************************/

#ifndef Geometryhpp
#define Geometryhpp

#include "Const.hpp"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cstring> 
using namespace std;

void ShowMyvector(myvector &cc)
{
   for(int i=0; i<cc.size(); ++i)
   {
      cout<<cc[i]<<'\t';
   }
   cout<<endl;
}

//Functions using vectors
inline void crossproduct(myvector &c1, myvector &c2, myvector &cc)
{
   cc.at(0) = c1.at(1) * c2.at(2) - c1.at(2) * c2.at(1);
   cc.at(1) = c1.at(2) * c2.at(0) - c1.at(0) * c2.at(2);
   cc.at(2) = c1.at(0) * c2.at(1) - c2.at(0) * c1.at(1);
}

inline float innerproduct(myvector &c1, myvector &c2)
{
   return c1.at(0) * c2.at(0) + c1.at(1) * c2.at(1) + c1.at(2) * c2.at(2);
}

inline void vectorsum(myvector &c1, myvector &c2, myvector &cc)
{
   cc.at(0) = c1.at(0) + c2.at(0);
   cc.at(1) = c1.at(1) + c2.at(1);
   cc.at(2) = c1.at(2) + c2.at(2);
}

inline void subtract(myvector &c1, myvector &c2, myvector &cc)
{
   cc.at(0) = c1.at(0) - c2.at(0);
   cc.at(1) = c1.at(1) - c2.at(1);
   cc.at(2) = c1.at(2) - c2.at(2);
}

inline bool norm(myvector &c, myvector &cc)
{
   float len = sqrt(c.at(0)*c.at(0) + c.at(1)*c.at(1) +c.at(2)*c.at(2));
   if(len<Extra) return false;
   cc.at(0) = c.at(0)/len;
   cc.at(1) = c.at(1)/len;
   cc.at(2) = c.at(2)/len;
   return true;
}

inline float vectorlength(myvector &c)
{
   return sqrt(c.at(0)*c.at(0)+ c.at(1)*c.at(1) + c.at(2)*c.at(2));
}

inline void multi(float coefficient, myvector &c, myvector &cc)
{
   cc.at(0) = c.at(0)*coefficient;
   cc.at(1) = c.at(1)*coefficient;
   cc.at(2) = c.at(2)*coefficient;
}

inline float deg2rad(float deg)
//float deg2rad(float deg)
{
   return deg*PI/180;
}

inline float rad2deg(float rad)
//float rad2deg(float rad)
{
   return rad*180/PI;
}


inline float Points2Distance2(myvector &c1, myvector &c2)
{
   return (c1.at(0)-c2.at(0))*(c1.at(0)-c2.at(0))
   +(c1.at(1)-c2.at(1))*(c1.at(1)-c2.at(1))
   +(c1.at(2)-c2.at(2))*(c1.at(2)-c2.at(2));
}

inline float Points2Distance(myvector &c1, myvector &c2)
{
   return sqrt(
   (c1.at(0)-c2.at(0))*(c1.at(0)-c2.at(0))
   +(c1.at(1)-c2.at(1))*(c1.at(1)-c2.at(1))
   +(c1.at(2)-c2.at(2))*(c1.at(2)-c2.at(2)) );
}

//Return the angle of c1-c2-c3. Unit: radian
//Angle <c1c2c3 
inline float Points2Angle(myvector &c1, myvector &c2, myvector &c3)
{
   float a=(c1.at(0)-c2.at(0))*(c1.at(0)-c2.at(0))
   +(c1.at(1)-c2.at(1))*(c1.at(1)-c2.at(1))
   +(c1.at(2)-c2.at(2))*(c1.at(2)-c2.at(2));
   
   float b=(c2.at(0)-c3.at(0))*(c2.at(0)-c3.at(0))
   +(c2.at(1)-c3.at(1))*(c2.at(1)-c3.at(1))
   +(c2.at(2)-c3.at(2))*(c2.at(2)-c3.at(2));
   
   float c=(c3.at(0)-c1.at(0))*(c3.at(0)-c1.at(0))
   +(c3.at(1)-c1.at(1))*(c3.at(1)-c1.at(1))
   +(c3.at(2)-c1.at(2))*(c3.at(2)-c1.at(2));
   
   float alpha=(a+b-c)/(2*sqrt(a*b));
   
   if (alpha>1 && alpha<1+Extra)
   {
      alpha=1;
   }
   else if(alpha<-1 && alpha>-1-Extra)
   {
      alpha=-1;
   }
   else if(alpha>1+Extra || alpha<-1-Extra)
   {
      cout<<"Error, float Points2Angle()\n";
      exit(0);
   }
   
   return acos(alpha);
}


/************************Points to Dihedral************************/
//If use inline here there will be a error when combiling codes.
//Return the angle of c1-c2-c3-c4. Unit: radian
//points c1-c2-c3-c4 should be in two different pane.
//or there may be faults.
inline float Points2Dihedral(myvector &c1, myvector &c2, myvector &c3, myvector &c4)
{
   myvector vector1(3,0), vector2(3,0), vector3(3,0);
   subtract(c1, c2, vector1);
   subtract(c2, c3, vector2);
   subtract(c3, c4, vector3);
   
   myvector v1(3,0), v2(3,0);
   crossproduct(vector2, vector1, v1);
   crossproduct(vector3, vector2, v2);
   
   myvector v3(3,0), v4(3,0);
   if(!norm (v1, v3)) {cout<<"Error in Points2Dihedral 1\n"<<endl; exit(1);}
   if(!norm (v2, v4)) {cout<<"Error in Points2Dihedral 2\n"<<endl; exit(1);}
   
   float dihedral = innerproduct(v3, v4);
   
   if (dihedral>1 && dihedral<1+Extra)
   {
      //cout<<"dihedral "<<dihedral<<" "<<acos(dihedral)<<"\n";
      dihedral=1;
   }
   else if(dihedral<-1 && dihedral>-1-Extra)
   {
      dihedral=-1;
   }
   else if(dihedral>1+Extra || dihedral<-1-Extra)
   {
      cout<<"Error, float Points2Dihedral()\n";
      exit(0);
   }
   
   myvector v5(3,0);
   crossproduct(v4, v3, v5);
   float direction = innerproduct(v5, vector2);
   
   if (direction>0)
    {
       return  acos(dihedral);
    }  
    else
    {  
       return -acos(dihedral);
    }
   
}

inline bool Points2Dihedral(myvector &c1, myvector &c2, myvector &c3, myvector &c4, float &dihedral)
{
   myvector vector1(3,0), vector2(3,0), vector3(3,0);
   subtract(c1, c2, vector1);
   subtract(c2, c3, vector2);
   subtract(c3, c4, vector3);
   
   myvector v1(3,0), v2(3,0);
   crossproduct(vector2, vector1, v1);
   crossproduct(vector3, vector2, v2);
   
   myvector v3(3,0), v4(3,0);
   if(!norm (v1, v3)) 
   {
      cout<<"Error in Points2Dihedral 1\n"<<endl; 
      return false;
   }
   if(!norm (v2, v4)) 
   {
      cout<<"Error in Points2Dihedral 2\n"<<endl; 
      return false;
   }
   
   dihedral = innerproduct(v3, v4);
   
   if (dihedral>1 && dihedral<1+Extra)
   {
      //cout<<"dihedral "<<dihedral<<" "<<acos(dihedral)<<"\n";
      dihedral=1;
   }
   else if(dihedral<-1 && dihedral>-1-Extra)
   {
      dihedral=-1;
   }
   else if(dihedral>1+Extra || dihedral<-1-Extra)
   {
      cout<<"Error, float Points2Dihedral()\n";
      exit(0);
   }
   
   myvector v5(3,0);
   crossproduct(v4, v3, v5);
   float direction = innerproduct(v5, vector2);
   
   if (direction>0)
    {
       dihedral=acos(dihedral);
       return true;
    }  
    else
    {  
       dihedral=-acos(dihedral);
       return true;
    }
   
}
//Given 2 points c1,c2 in the 3D space, generate the third points on the 
//line of c1, c2, with a distance of dis32 to c2;
//Date 2008.5.28.
//BUG: PointsOnLine()
//c3[0]=c2[0]-(c1[0]-c2[0])*dis32/dis12; => c3[0]=c2[0]-(c1[0]-c2[0])*dis32/dis12;
//Date 2008.7.11.
inline bool PointsOnLine(float dis32, myvector &c1, myvector &c2, myvector &c3)
{
   float dis12=Points2Distance(c1, c2);
   if(dis12<Extra) return false;
   
   c3.assign(3,0);
   c3[0]=c2[0]-(c1[0]-c2[0])*dis32/dis12;
   c3[1]=c2[1]-(c1[1]-c2[1])*dis32/dis12;
   c3[2]=c2[2]-(c1[2]-c2[2])*dis32/dis12;
   return true;
}

/******************************************************************************************************************
                          Function: internal2cartesian
                                       by 
                     CY, Jiang' Lab. Institute of Biophysics
                                             2005.10.18.

Calculate the Cartesian coordinates from internal coordinates.

Parameter : 3 Points' Cartesian coordinates (c1[3], c2[3], c3[3]) and 1 point's internal 
            coordinats ( p[3] ). Parameter c4[3] is p[3]'s cartesian coordinates calculated
            by the function.
            p[0]: distance, p[1]: angle/deg, p[2]=dihedral/deg

Return bool variable.
True  : successful calculation.
False : fails. The reason is that the 3 points are in one line. 

******************************************************************************************************************/

bool internal2cartesian(myvector &c1, myvector &c2, myvector &c3, myvector &p, myvector &c4)
{
   myvector d1(3,0), d2(3,0) ,xp(3,0);

   subtract(c1, c2, d1);
   subtract(c3, c2, d2);
   crossproduct(d1, d2, xp);

   if( (xp[0]<Extra && xp[0]>-Extra) && (xp[1]<Extra && xp[1]>-Extra)
    && (xp[2]<Extra && xp[2]>-Extra) )
    {
        cout<<"Error! Points 1, 2, 3 are in one line & no plane!\n";
        ShowMyvector(c1);ShowMyvector(c2);ShowMyvector(c3);
        return false;        
    }

   myvector d3(3,0), yp(3,0), r(3,0), ypp(3,0), tmp1(3,0), tmp2(3,0);

   crossproduct(d2, xp, yp);
   float ang1 = deg2rad(p.at(2));
   norm(xp, xp);
   norm(yp, yp);
   multi(cos(ang1), xp, tmp1);
   multi(sin(ang1), yp, tmp2);
   vectorsum(tmp1, tmp2, r);

   crossproduct(d2, r, ypp);
   float ang2 = deg2rad(p.at(1));
   norm(d2, d2 );
   norm(ypp,ypp);
   multi(-cos(ang2), d2, tmp1);
   multi(sin(ang2), ypp, tmp2);
   vectorsum(tmp1, tmp2, d3);
   multi(p.at(0), d3, d3);
   vectorsum(c3, d3, d3);

   c4.assign(3, 0);
   c4.at(0) = d3[0];
   c4.at(1) = d3[1];
   c4.at(2) = d3[2];

   return true;
}





//Translate cartesian to internal coordinate system. 
//Unit: degree
//p[0]: distance, p[1]: angle, p[2]=dihedral
bool cartesian2internal (myvector &c1, myvector &c2, myvector &c3, myvector &c4, myvector &p)
{
   if(c1.size()!=3 || c2.size()!=3 || c3.size()!=3 || c4.size()!=3)
      return false;
   p.assign(3, 0);
   p[0]=Points2Distance(c3, c4);
   p[1]=rad2deg( Points2Angle(c2, c3, c4) );
   p[2]=rad2deg( Points2Dihedral(c1, c2, c3, c4) );
   
   return true;
}






//Functions using arrays

inline void crossproduct(float *c1, float *c2, float *cc)
{
   cc[0] = c1[1] * c2[2] - c1[2] * c2[1];
   cc[1] = c1[2] * c2[0] - c1[0] * c2[2];
   cc[2] = c1[0] * c2[1] - c2[0] * c1[1];
}

inline float innerproduct(float *c1, float *c2)
{
   return c1[0] * c2[0] + c1[1] * c2[1] + c1[2] * c2[2];
}

inline void vectorsum(float *c1, float *c2, float *cc)
{
   cc[0] = c1[0] + c2[0];
   cc[1] = c1[1] + c2[1];
   cc[2] = c1[2] + c2[2];
}

inline void subtract(float *c1, float *c2, float *cc)
{
   cc[0] = c1[0] - c2[0];
   cc[1] = c1[1] - c2[1];
   cc[2] = c1[2] - c2[2];
}

inline bool norm(float* c, float* cc)
{
   float len = sqrt(c[0]*c[0] + c[1]*c[1] +c[2]*c[2]);
   if(len<Extra) return false;
   cc[0] = c[0]/len;
   cc[1] = c[1]/len;
   cc[2] = c[2]/len;
   return true;
}

inline void multi(float coefficient, float *c)
{
   c[0] *= coefficient;
   c[1] *= coefficient;
   c[2] *= coefficient;
}

/************************Points to Dihedral************************/
//If use inline here there will be a error when combiling codes.
float Points2Dihedral(float *c1, float *c2, float *c3, float *c4)
{
   float vector1[3], vector2[3], vector3[3];
   float v1[3], v2[3], v3[3], v4[3], v5[3];
   
   subtract(c1, c2, vector1);
   subtract(c2, c3, vector2);
   subtract(c3, c4, vector3);
   
   crossproduct(vector2, vector1, v1);
   crossproduct(vector3, vector2, v2);
   
   norm (v1, v3);
   norm (v2, v4);
   
   float dihedral = innerproduct(v3, v4);
   
   if (dihedral>1 && dihedral<1+Extra)
   {
      //cout<<"dihedral "<<dihedral<<" "<<acos(dihedral)<<"\n";
      dihedral=1;
   }
   else if(dihedral<-1 && dihedral>-1-Extra)
   {
      dihedral=-1;
   }
   else if(dihedral>1+Extra || dihedral<-1-Extra)
   {
      cout<<"Error, float Points2Dihedral()\n";
      exit(0);
   }
   
   crossproduct(v4, v3, v5);
   float direction = innerproduct(v5, vector2);
   
   if (direction>0)
    {
       return  acos(dihedral);
    }  
    else
    {  
       return -acos(dihedral);
    }
   
}



/***************************************************************
//Translate Internal data into Cartesian data for 
//	the last of the four points
***************************************************************/
bool internal2cartesian (float *c1, float *c2, float *c3, float *p, float *c4)
{
   float d1[3];
   float d2[3];
   float xp[3];

   subtract(c1, c2, d1);
   subtract(c3, c2, d2);
   crossproduct(d1, d2, xp);

   if( (xp[0]<Extra && xp[0]<-Extra) && (xp[1]<Extra && xp[1]>-Extra) && (xp[2]<Extra && xp[2]>-Extra) )
    {
        cout<<"Error! Points 1, 2, 3 are in one line & no plane!\n";
        return false;        
    }

   float d3[3];
   float yp[3];
   float r[3] ;
   float ypp[3];
   float tmp1[3];

   crossproduct(d2, xp, yp);
   float ang1 = deg2rad(p[2]);
   norm(xp, xp);
   norm(yp, xp);
   multi(cos(ang1), xp);
   multi(sin(ang1), yp);
   vectorsum( xp, yp, r);

   crossproduct(d2, r, ypp);
   float ang2 = deg2rad(p[1]);
   norm(d2, d2 );
   norm(ypp, ypp);
   multi(-cos(ang2), d2);
   multi(sin(ang2), ypp);
   vectorsum( d2, ypp, d3 );
   multi(p[0], d3);
   vectorsum(c3, d3, tmp1);

   c4[0] = tmp1[0];
   c4[1] = tmp1[1];
   c4[2] = tmp1[2];

   return true;
}


#endif
