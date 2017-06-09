//能量计算部分从OptSite18.cpp抽提出来，其他来自SolHbond6.cpp
//改来可以计算单独残基的溶剂能2009.3.5.
//g++ wesol2.cpp -o wesol2
//caoyang


# include "SimpleProteinTools.hpp"
# include "Rotamer2Structure4.hpp"
# include "simple_water_model.hpp"
# include "StericOccupy3.hpp"

using namespace std;

float WaterRadii=1.4;
float latticesize=0.8; //格点空间单位长度
float _cara=5; // 计算表面自由空间半径
float enlar=1; // 表面积权重
float hbsphere=3;  //表面氢键有效半径

float Kentropy=1.0;
float Khbonds1=-2.0;
float Khbonds2=-1.5;

char AAname3to1(string name3)
{
	char name1;
	if( name3=="ALA" )	name1='A';
	else if( name3=="ARG" ) name1='R';
	else if( name3=="ASN" ) name1='N';
	else if( (name3=="ASP")||(name3=="ASX") ) name1='D';
	else if( name3=="CYS" ) name1='C';
	else if( (name3=="GLN")||(name3=="GLX") ) name1='Q';
	else if( name3=="GLU" ) name1='E';
	else if( name3=="GLY" ) name1='G';
	else if( name3=="HSD" ) name1='H';
	else if( (name3=="HIS")||(name3=="HSE")||(name3=="HSD")||(name3=="HSP") ) name1='H';
	else if( name3=="ILE" ) name1='I';
	else if( name3=="LEU" ) name1='L';
	else if( name3=="LYS" ) name1='K';
	else if( name3=="MET" ) name1='M';
	else if( name3=="PHE" ) name1='F';
	else if( name3=="PRO" ) name1='P';
	else if( name3=="SER" ) name1='S';
	else if( name3=="THR" ) name1='T';
	else if( name3=="TRP" ) name1='W';
	else if( name3=="TYR" ) name1='Y';
	else if( name3=="VAL" ) name1='V';
	else
	{
		name1='*';
		//cout << "\nAAname3to1()--Error: no matched residue name.\t" << name3 << endl;
	}
	
	return name1;
}

   
static float AtomRadius(char Name)
{
   if(Name=='C') return 2.0; //1.7; 
   if(Name=='N') return 1.85;//1.6; 
   if(Name=='O') return 1.7; //1.4; 
   if(Name=='S') return 2.0; //1.8; 
   if(Name=='H') return 1.2; 
   else return 1.0;
}
float SolEn(strs &peptide, array2sd &cartesian, vector<float> &hbes)
{  
   static int i=0, j=0, k=0, l=0;
   static float radii=WaterRadii;
   float square_unit=latticesize*latticesize;
   
   static str2vectord::iterator pos2;
   vector<vector<float> > points, bound, elepoi;
   vector<float> elepot;
   vector<int> atom2residue;
   int atmele=0;
   for(int i=0; i<peptide.size(); i++)
   {
      for(pos2=cartesian[i].begin(); pos2!=cartesian[i].end(); ++pos2)
      {
         l=SphereBuilding(pos2->second, AtomRadius(pos2->first[0])+radii, latticesize, points);
         for(k=0; k<l; k++)
            atom2residue.push_back(i);
      }
   }
   StericOccupy so(latticesize, points);
   //so.ShowSize();

   vector<int> effpoint, effbound, bound2residue;
   vector<vector<int> > atomindex;
   int amount=so.VertexAmount(effpoint, atomindex);
   //cout<<"Amount "<<amount<<endl;
   for(i=0; i<effpoint.size(); i++)
   {
      if(so.BoundaryFaceIndex(effpoint[i]))
      {
         effbound.push_back(effpoint[i]);
         bound2residue.push_back(atom2residue.at(atomindex[i][0]));
         //for(k=0; k<atomindex[i].size(); k++)
         //cout<<peptide[atom2residue.at(atomindex[i][k])]<<" ";
         //cout<<"\n";
      }
   }
   float tempfree=0, logspace=0, sumspace=0, sdensity=0, sumdense=0, spotenti=0, sumpoten=0;
   vector<int> max3;
   so.CopyLattice(max3);
   StericOccupy effsu(latticesize, max3, effbound);
   //effsu.ShowSize();
   
   vector<float> reslogspace(peptide.size(), 0);
   vector<float> reshydrogen(peptide.size(), 0);
   
   for(i=0; i<effbound.size(); i++)
   {
      tempfree=so.FreeSpaceRatio(effbound[i], _cara);
      sdensity=1/((float)effsu.NeighborIndex(effbound[i], (int)1));
      tempfree=-log(tempfree)*sdensity;
      logspace+=tempfree;
      reslogspace[bound2residue[i]]+=tempfree*square_unit;
   }
   
   SetWR(radii);
   static vector<myvector> wats; 
   static vector<int> param; 
   vector<float> reshb(7, 0);
   hbes.assign(7, 0);
   
   vector<int> num_waters;
   
   for(i=0; i<peptide.size(); i++)
   {
      wats.clear(); param.clear(); reshb.assign(7, 0);
      PlaceWater(peptide[i], cartesian[i], wats, param);
      for(j=0; j<wats.size(); j++)
      {
         tempfree=so.FreeSpaceRatio(wats[j], hbsphere);
         hbes[param[j]-1]+=tempfree;
         reshb.at(param[j]-1)+=tempfree;
      }
      num_waters.push_back(wats.size());
      reshydrogen[i]= reshb[0]*Khbonds1+reshb[1]*Khbonds2+reshb[2]*Khbonds2
                     +reshb[3]*Khbonds1+reshb[4]*Khbonds2+reshb[5]*Khbonds2
                     +reshb[6]*Khbonds2;
   }
   for(i=0; i<peptide.size(); i++)
   {
      //cout<<i<<" "<<AAname3to1(peptide[i])<<" "<<reslogspace[i]+reshydrogen[i]<<endl;
      cout<<i<<" "<<AAname3to1(peptide[i])<<" "<<reshydrogen[i] << ' ' << reslogspace[i] <<endl;
   }
   
   
   
   return logspace*square_unit;
}

static float SolvationFreeEnergy(strs &peptide, array2sd &temp_coor)
{
      static vector<float> hbs;
      float ens=SolEn(peptide, temp_coor, hbs);
      //cout<<"Ens ";
      //cout<<ens<<" ";
      //for(int i=0; i<7; i++)
      {
         //cout<<hbs[i]<<" ";
      }
      //cout<<"\n";
      
      float solvation_energy=0;
      solvation_energy+=ens;
      
      solvation_energy+=hbs[0]*Khbonds1+hbs[1]*Khbonds2+hbs[2]*Khbonds2
                       +hbs[3]*Khbonds1+hbs[4]*Khbonds2+hbs[5]*Khbonds2
                       +hbs[6]*Khbonds2;
      
      return solvation_energy;
}



int main(int argc, char **argv)
{
   if(argc<2)
   {
   cout<<"Error input!\nUsage: wesol2 [Input PDB] [格点模型单位边长]\n"; 
   cout<<"e.g. ./wesol2 1R69.pdb 0.7\n";

   exit(0);
   }
   char *input1=*(argv+1);
   ProteinStructureAnalysis psa1;
   psa1.GetData(input1);
   strs peptide; array2sd cartesian;
   psa1.ExtractPDB(&peptide, &cartesian);
   psa1.RemoveHydrogen(cartesian);
   
   if(argc>2)  latticesize=atof(argv[2]);
   if(argc>3)  WaterRadii=atof(argv[3]);
   if(argc>4)  hbsphere=atof(argv[4]);
   if(argc>5)  Kentropy=atof(argv[5]);
   if(argc>6)  Khbonds1=atof(argv[6]);
   if(argc>7)  Khbonds2=atof(argv[7]);
   
   Rotamer2Structure r2s(top_para.c_str());
      for(int i=0, k=0; i<peptide.size(); i++)
      {
         if(peptide[i].size()==0) {k=i+1; continue;}
         if(i>0 && i!=k)  r2s.ConstructHN(peptide[i], cartesian[i-1], cartesian[i]);
         r2s.ConstructHs(peptide[i], cartesian[i]); 
      }
   //cout<<"Build Hydrogen Finished...\n";
   //cout<<"Solv: "<<SolvationFreeEnergy(peptide, cartesian)<<endl;
   SolvationFreeEnergy(peptide, cartesian);
   return 0;
}



