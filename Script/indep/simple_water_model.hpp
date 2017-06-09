# include "Geometry.hpp"
//Author caoyang
using namespace std;

float wr=1.4;
float SetWR(float swr)
{
   wr=swr;
   return wr;
}

//Carbonyl Oxygen
//c1: O; c2: C; C2: other atom
int PositionWater1(myvector &c1, myvector &c2, myvector &c3, vector<myvector> &ps)
{
   int i=0, j=0, k=0, n=0;
   myvector _internal(3, 0), temp, h0, h1, h2 ;
   h0.push_back(1.7+wr);
   h1.push_back(120);
   h2.push_back(0);
   h2.push_back(180);
   
   for(i=0; i<h2.size(); i++)
   {
      _internal[2]=h2[i];
      for(j=0; j<h1.size(); j++)
      {
         _internal[1]=h1[j];
         for(k=0; k<h0.size(); k++)
         {
            _internal[0]=h0[k];
            if(!internal2cartesian(c3, c2, c1, _internal, temp))
            {
               cout<<"Error in PositionWater1\n";
               break;
            }
            ps.push_back(temp);
            n++;
         }
      }
   }
   return n;
}

//Amide nitrogen
//c1: H, c2: N
int PositionWater2(myvector &c1, myvector &c2, vector<myvector> &ps)
{
   int k=0, n=0;
   myvector temp(3, 0), h0;
   //h0.push_back(1.85+wr);
   h0.push_back(1.2+wr);
   
   for(k=0; k<h0.size(); k++)
   {
      if(!PointsOnLine(h0[k], c2, c1, temp))
      {
         cout<<"Error in PositionWater2\n";
         break;
      }
      ps.push_back(temp);
      n++;
   }
   return n;
}
//Hydroxyl oxygen
//For OH in SER THR: (c1, c2, c3)顺序是：H-C-O， H-C不成键
int PositionWater3(myvector &c1, myvector &c2, myvector &c3, vector<myvector> &ps)
{
   int i=0, j=0, k=0, n=0;
   myvector _internal(3, 0), temp(3, 0), h0, h1, h2 ;
   h0.push_back(1.7+wr);
   h1.push_back(110);
   h2.push_back(120);
   h2.push_back(240);

   for(i=0; i<h2.size(); i++)
   {
      _internal[2]=h2[i];
      for(j=0; j<h1.size(); j++)
      {
         _internal[1]=h1[j];
         for(k=0; k<h0.size(); k++)
         {
            _internal[0]=h0[k];
            if(!internal2cartesian(c1, c2, c3, _internal, temp))
            {
               cout<<"Error in PositionWater3 A\n";
               break;
            }
            ps.push_back(temp);
            n++;
         }
      }
   }
   h0.clear();
   h0.push_back(1.2+wr);
   for(k=0; k<h0.size(); k++)
   {
      if(!PointsOnLine(h0[k], c3, c1, temp))
      {
         cout<<"Error in PositionWater3 B\n";
         break;
      }
      ps.push_back(temp);
      n++;
   }
   
   return n;
}
//Hydroxyl oxygen
//For TYR (c1, c2, c3)顺序是：H-C-O， H-C不成键
int PositionWater4(myvector &c1, myvector &c2, myvector &c3, vector<myvector> &ps)
{
   int i=0, j=0, k=0, n=0;
   myvector _internal(3, 0), temp(3, 0), h0, h1, h2 ;
   h0.push_back(1.7+wr);
   h1.push_back(110);
   h2.push_back(180);
   
   for(i=0; i<h2.size(); i++)
   {
      _internal[2]=h2[i];
      for(j=0; j<h1.size(); j++)
      {
         _internal[1]=h1[j];
         for(k=0; k<h0.size(); k++)
         {
            _internal[0]=h0[k];
            if(!internal2cartesian(c1, c2, c3, _internal, temp))
            {
               cout<<"Error in PositionWater4 A\n";
               break;
            }
            ps.push_back(temp);
            n++;
         }
      }
   }
   h0.clear();
   h0.push_back(1.2+wr);
   for(k=0; k<h0.size(); k++)
   {
      if(!PointsOnLine(h0[k], c3, c1, temp))
      {
         cout<<"Error in PositionWater4 B\n";
         break;
      }
      ps.push_back(temp);
      n++;
   }
   return n;
}

//Aromatic nitrogen
//For HIS (c1, c2, c3)顺序是：N-N-H， N-N不成键
int PositionWater5(myvector &c1, myvector &c2, myvector &c3, vector<myvector> &ps)
{
   int k=0, n=0;
   myvector _internal(3, 0), temp(3, 0), h0;
   
   h0.push_back(1.85+wr);
   for(k=0; k<h0.size(); k++)
   {
      if(!PointsOnLine(h0[k], c2, c1, temp))
      {
         cout<<"Error in PositionWater5 A\n";
         break;
      }
      ps.push_back(temp);
      n++;
   }
   h0.clear();
   h0.push_back(1.2+wr);
   for(k=0; k<h0.size(); k++)
   {
      if(!PointsOnLine(h0[k], c2, c3, temp))
      {
         cout<<"Error in PositionWater5 B\n";
         break;
      }
      ps.push_back(temp);
      n++;
   }
   
   
   return n;
}
/*
int WaterType(const char type, myvector & param)
{
   switch (type)
   //asp, glu O
   case 'A': param.push_back(2.0); return 1;
   //=O of backbone and asn, gln
   case 'B': param.push_back(1.5); return 2;
   //O of hydroxyl
   case 'C': param.push_back(1.0); return 3;
   //HN of amid nitrogen
   case 'D': param.push_back(2.0); return 4;
   //N of HIS
   case 'E': param.push_back(1.5); return 5;
   default: cout<<"Error water type\n"; return 0;
}
*/


int PlaceWater(string &res, str2vectord &coor, vector<myvector> &wats, vector<int> &param)
{
   if(coor.size()==0) return 0;
   int num=0;
   if(res=="ASP")
   {
      num+=PositionWater1(coor["OD1"], coor["CG"], coor["CB"], wats);
      num+=PositionWater1(coor["OD2"], coor["CG"], coor["CB"], wats);
      
      param.push_back(1); param.push_back(1); 
      param.push_back(1); param.push_back(1); 
   }
   else if(res=="ASN")
   {
      num+=PositionWater1(coor["OD1"], coor["CG"], coor["CB"], wats);
      num+=PositionWater2(coor["HD21"], coor["ND2"], wats);
      num+=PositionWater2(coor["HD22"], coor["ND2"], wats);
      
      param.push_back(2); param.push_back(2); 
      param.push_back(3); param.push_back(3); 
   }
   else if(res=="ARG")
   {
      num+=PositionWater2(coor["HH11"], coor["NH1"], wats);
      num+=PositionWater2(coor["HH12"], coor["NH1"], wats);
      num+=PositionWater2(coor["HH21"], coor["NH2"], wats);
      num+=PositionWater2(coor["HH22"], coor["NH2"], wats);
      num+=PositionWater2(coor["HE"], coor["NE"], wats);
      
      param.push_back(4); param.push_back(4);
      param.push_back(4); param.push_back(4);
      param.push_back(3);
   }
   else if(res=="LYS")
   {
      num+=PositionWater2(coor["HZ1"], coor["NZ"], wats);
      num+=PositionWater2(coor["HZ2"], coor["NZ"], wats);
      num+=PositionWater2(coor["HZ3"], coor["NZ"], wats);
      
      param.push_back(4);param.push_back(4);param.push_back(4);
   }
   else if(res=="GLU")
   {
      num+=PositionWater1(coor["OE1"], coor["CD"], coor["CG"], wats);
      num+=PositionWater1(coor["OE2"], coor["CD"], coor["CG"], wats);
      
      param.push_back(1); param.push_back(1); 
      param.push_back(1); param.push_back(1); 
   }
   else if(res=="GLN")
   {
      num+=PositionWater1(coor["OE1"], coor["CD"], coor["CG"], wats);
      num+=PositionWater2(coor["HE21"], coor["NE2"], wats);
      num+=PositionWater2(coor["HE22"], coor["NE2"], wats);
      
      param.push_back(2); param.push_back(2); 
      param.push_back(3); param.push_back(3); 
   }
   else if(res=="SER")
   {
      //For OH in SER THR: (c1, c2, c3)顺序是：H-C-O， H-C不成键
      num+=PositionWater3(coor["HG1"], coor["CB"], coor["OG"], wats);
      
      param.push_back(5); param.push_back(5);
      param.push_back(6);
   }
   else if(res=="THR")
   {
      num+=PositionWater3(coor["HG1"], coor["CB"], coor["OG1"], wats);
      
      param.push_back(5); param.push_back(5);
      param.push_back(6);
   }
   else if(res=="TYR")
   {
      //For TYR (c1, c2, c3)顺序是：H-C-O， H-C不成键
      num+=PositionWater4(coor["HH"], coor["CZ"], coor["OH"], wats);
      
      param.push_back(5); param.push_back(6);
   }
   else if(res=="TRP")
   {
      num+=PositionWater2(coor["HE1"], coor["NE1"], wats);
      
      param.push_back(3);
   }
   else if(res=="HSD")
   {
      //For HIS (c1, c2, c3)顺序是：N-N-H， N-N不成键
      num+=PositionWater5(coor["NE2"], coor["ND1"], coor["HD1"], wats);
      
      param.push_back(7); param.push_back(3);
   }
   else if(res=="SAB" || res=="SAC")
   {
      num+=PositionWater1(coor["O1A"], coor["C1"], coor["C2"], wats);
      num+=PositionWater1(coor["O1B"], coor["C1"], coor["C2"], wats);
      num+=PositionWater3(coor["HO4"], coor["C4"], coor["O4"], wats);
      num+=PositionWater2(coor["HN"], coor["N5"], wats);
      num+=PositionWater1(coor["O10"], coor["C10"], coor["C11"], wats);
      num+=PositionWater3(coor["HO7"], coor["C7"], coor["O7"], wats);
      num+=PositionWater3(coor["HO8"], coor["C8"], coor["O8"], wats);
      num+=PositionWater3(coor["HO9"], coor["C9"], coor["O9"], wats);
      
      param.push_back(2); param.push_back(2); 
      param.push_back(2); param.push_back(2);
      param.push_back(5); param.push_back(5); param.push_back(6);
      param.push_back(3);
      param.push_back(2); param.push_back(2); 
      param.push_back(5); param.push_back(5); param.push_back(6);
      param.push_back(5); param.push_back(5); param.push_back(6);
      param.push_back(5); param.push_back(5); param.push_back(6);
   }
   else if(res=="GLB")
   {
      num+=PositionWater3(coor["HO2"], coor["C2"], coor["O2"], wats);
      num+=PositionWater3(coor["HO3"], coor["C3"], coor["O3"], wats);
      num+=PositionWater3(coor["HO4"], coor["C4"], coor["O4"], wats);
      
      param.push_back(5); param.push_back(5); param.push_back(6);
      param.push_back(5); param.push_back(5); param.push_back(6);
      param.push_back(5); param.push_back(5); param.push_back(6);
   }
   else if(res=="GLD")
   {
      num+=PositionWater3(coor["HO1"], coor["C1"], coor["O1"], wats);
      num+=PositionWater3(coor["HO2"], coor["C2"], coor["O2"], wats);
      num+=PositionWater3(coor["HO3"], coor["C3"], coor["O3"], wats);
      num+=PositionWater3(coor["HO4"], coor["C4"], coor["O4"], wats);
      
      param.push_back(5); param.push_back(5); param.push_back(6);
      param.push_back(5); param.push_back(5); param.push_back(6);
      param.push_back(5); param.push_back(5); param.push_back(6);
      param.push_back(5); param.push_back(5); param.push_back(6);
   }
   else if(res=="GLC")
   {
      num+=PositionWater3(coor["HO2"], coor["C2"], coor["O2"], wats);
      num+=PositionWater3(coor["HO4"], coor["C4"], coor["O4"], wats);
      num+=PositionWater3(coor["HO6"], coor["C6"], coor["O6"], wats);
      
      param.push_back(5); param.push_back(5); param.push_back(6);
      param.push_back(5); param.push_back(5); param.push_back(6);
      param.push_back(5); param.push_back(5); param.push_back(6);
   }
   else if(res=="NGB" || res=="NGC")
   {
      num+=PositionWater2(coor["HN"], coor["N2"], wats);
      num+=PositionWater3(coor["HO3"], coor["C3"], coor["O3"], wats);
      num+=PositionWater3(coor["HO6"], coor["C6"], coor["O6"], wats);
      num+=PositionWater1(coor["O7"], coor["C7"], coor["C8"], wats);
      
      param.push_back(3);
      param.push_back(5); param.push_back(5); param.push_back(6);
      param.push_back(5); param.push_back(5); param.push_back(6);
      param.push_back(2); param.push_back(2); 
   }
   
   if(coor.count("HN")==1 && coor.count("N")==1)
   {
      num+=PositionWater2(coor["HN"], coor["N"], wats);
      
      param.push_back(3);
   }
   if(coor.count("O")==1 && coor.count("C")==1 && coor.count("CA")==1 )
   {
      num+=PositionWater1(coor["O"], coor["C"], coor["CA"], wats);
      
      param.push_back(2); param.push_back(2); 
   }
   //cout<<"Generate Water Number "<<num<<endl;
   return num;
}

