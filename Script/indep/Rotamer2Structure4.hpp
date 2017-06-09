//Build side chain from CB
//The old one is from CA
//Date 2007.11.2.
//Both CA, CB is OK. Just using define CB in Const.hpp or not. 
//Date 2007.11.21.
//Add int ConstructHs(const string& res, str2vectord &group); 
//and Modified OutPdbFile() into pure PDB format
//Date 2008.4.25.
//void AtomIndexMaping => int AtomIndexMaping
//int ConstructHN and int ConstructHs return 0 if there is no query res type
//Date 2008.5.8.
//int OutPdbFileALL() for any kind of atoms
//Date 2008.7.10
//chain id is null and start number is 1
//Add HSP, so the residue type is 21. which is different from version3
//Date 2008.11.3.
#ifndef Rotamer2Structure2Head
#define Rotamer2Structure2Head

const int residue_type_num=21;

//#define discal // return distance between CB and the fartest side chain atom.

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <deque>
#include <iomanip>
#include <cstdlib>
#include <cstring> 
#include "Const.hpp"
#include "Geometry.hpp"

using namespace std;

int PDB_Format_OUT(int j, string atom, string residue, int i, vector<float> &coordinates)
{          
      cout<<setiosflags(ios::left)<<setw(6)<<"ATOM"<<resetiosflags(ios::left);
      cout<<setw(5)<<j<<" ";
      if(atom.size()==4)
      cout<<setiosflags(ios::left)<<setw(3)<<atom<<resetiosflags(ios::left);
      else
      cout<<" "<<setiosflags(ios::left)<<setw(3)<<atom<<resetiosflags(ios::left);
      cout<<" ";
      cout<<residue;
      cout<<" ";
      cout<<" "<<setw(4)<<i;
      cout<<"    ";
      cout<<setiosflags(ios::fixed);
      cout<<setprecision(3);
      if(coordinates.size()!=3) 
      {cout<<endl; return 1;}
      cout<<setw(8)<<coordinates.at(0);
      cout<<setw(8)<<coordinates.at(1);
      cout<<setw(8)<<coordinates.at(2);
      cout<<endl;
      
      return 1;
}

int PDB_Format_OUT(int j, string atom, string residue, int i, vector<float> &coordinates, ofstream &pdbout)
{          
      pdbout<<setiosflags(ios::left)<<setw(6)<<"ATOM"<<resetiosflags(ios::left);
      pdbout<<setw(5)<<j<<" ";
      if(atom.size()==4)
      pdbout<<setiosflags(ios::left)<<setw(3)<<atom<<resetiosflags(ios::left);
      else
      pdbout<<" "<<setiosflags(ios::left)<<setw(3)<<atom<<resetiosflags(ios::left);
      pdbout<<" ";
      pdbout<<residue;
      pdbout<<" ";
      pdbout<<" "<<setw(4)<<i;
      pdbout<<"    ";
      pdbout<<setiosflags(ios::fixed);
      pdbout<<setprecision(3);
      if(coordinates.size()!=3) 
      {pdbout<<endl; return 1;}
      pdbout<<setw(8)<<coordinates.at(0);
      pdbout<<setw(8)<<coordinates.at(1);
      pdbout<<setw(8)<<coordinates.at(2);
      pdbout<<endl;
      
      return 1;
}

int PDB_Format_OUT(int j, string atom, string residue, char chn, int i, vector<float> &coordinates, ofstream &pdbout)
{          
      pdbout<<setiosflags(ios::left)<<setw(6)<<"ATOM"<<resetiosflags(ios::left);
      pdbout<<setw(5)<<j<<" ";
      if(atom.size()==4)
      pdbout<<setiosflags(ios::left)<<setw(3)<<atom<<resetiosflags(ios::left);
      else
      pdbout<<" "<<setiosflags(ios::left)<<setw(3)<<atom<<resetiosflags(ios::left);
      pdbout<<" ";
      pdbout<<residue;
      pdbout<<" ";
      pdbout<<chn<<setw(4)<<i;
      pdbout<<"    ";
      pdbout<<setiosflags(ios::fixed);
      pdbout<<setprecision(3);
      if(coordinates.size()!=3) 
      {pdbout<<endl; return 1;}
      pdbout<<setw(8)<<coordinates.at(0);
      pdbout<<setw(8)<<coordinates.at(1);
      pdbout<<setw(8)<<coordinates.at(2);
      pdbout<<endl;
      
      return 1;
}

class Rotamer2Structure
{
public:
      Rotamer2Structure(const char* inp);
      int AtomIndexMaping(const string& res, str2vectord &coordata);
      
      float ConstructAASidechain(const string& res, str2vectord &NCACB, 
             const vector <float> &torsion, str2vectord &cartesian);
      float res2radius(string res);
      
      int Atom2Output(str2strs &outseq);
      int OutPdbFile(strs &peptide, array2sd &cartesian);
      int OutPdbFile(strs &peptide, array2sd &cartesian, ofstream &pdbout);
      int OutPdbFileALL(strs &peptide, array2sd &cartesian, ofstream &pdbout);
      int OutPdbFile(strs &peptide, array2sd &cartesian, vector<char> &chn, ofstream &pdbout);
      int OutPdbFileALL(strs &peptide, array2sd &cartesian, vector<char> &chn, ofstream &pdbout);
      
      int Switch2AmberParameter(const char* inp);
      //Add HN of main chain
      int ConstructHN(const string& res, str2vectord &mC, str2vectord &NCA);
      //Add H of side chain other than HN
      int ConstructHs(const string& res, str2vectord &group); 

private:
      strs index2res;
      str2int res2index;
      str2float res2radii;
      //atom to index of each residue
      str2int atom2index[residue_type_num];
      //default internal parameters
      matrix mainchainpara;
      matrix sidechainpara;
      matrix hydrogenspara;
      //default atom index in each residue for fast calculation
      vector<vector<int> > mainatomindex;
      vector<vector<int> > sideatomindex;
      vector<vector<int> > hydrogenindex;
      //default atom names of each residue
      vector<strs> mainatom;
      vector<strs> sideatom;
      vector<strs> hydrogen;
      
      matrix cycle;
      myvector xyz, intercoor, cartcoor;
      bool amber_para;
};
int Rotamer2Structure::Switch2AmberParameter(const char* inp)
{
      vector<int> index; strs name;
      amber_para=true;
      index2res.clear(); 
      res2index.clear(); 
      res2radii.clear(); 
      mainchainpara.assign(residue_type_num, xyz);
      sidechainpara.assign(residue_type_num, xyz);
      mainatomindex.assign(residue_type_num, index);
      sideatomindex.assign(residue_type_num, index);
      mainatom.assign(residue_type_num, name);
      sideatom.assign(residue_type_num, name);
      intercoor.assign(3, 0);
      cartcoor.assign(3, 0);
      
      ifstream oi(inp, ios::in);
      if(!oi.good()) cout<<"Error in opening para file: "<<inp<<endl;
      string res, mark;
      strs atoms(4);
      str2int temp;
      vector<float> data(5);
      int i=0, j=0, n=0;
      float radii=0;
      while(oi.good())
      {
         oi>>mark;
         if(mark=="RESI")
         {
            oi>>res>>radii;
            res2index[res]=i;
            res2radii[res]=radii;
            index2res.push_back(res);
            //cout<<i<<"   "<<res<<endl;
            j=3;
            temp.clear();
            temp["-C"]=0; mainatom[i].push_back("-C");
            temp["N"]=1; mainatom[i].push_back("CA");
            temp["CA"]=2;  mainatom[i].push_back("N");
         }
#ifdef CB
         else if(mark=="IC" || mark=="CC")
#else
         else if(mark=="IC")
#endif
         {
            oi>>atoms[0]>>atoms[1]>>atoms[2]>>atoms[3]>>data[0]>>data[1]>>data[2]>>data[3]>>data[4];
            if(atoms[2][0]=='*')
            {
               atoms[2].erase(0,1);
            }
            for(n=0; n<3; n++)
            {
               if(temp.count(atoms[n])==0)
               {
                  cout<<"Error 11amber "<<atoms[n]<<" "<<res<<endl;
                  exit(0);
               }
               mainatomindex[i].push_back(temp[atoms[n]]);
               mainchainpara[i].push_back(data[2+n]);
            }
            mainatom[i].push_back(atoms[3]);
            temp[atoms[3]]=j;
            j++;
         }
#ifdef CB
         else if(mark=="AB")
#else
         else if(mark=="AB" || mark=="CC")
#endif
         {
            oi>>atoms[0]>>atoms[1]>>atoms[2]>>atoms[3]>>data[0]>>data[1]>>data[2];
            if(atoms[2][0]=='*')
            {
               atoms[2].erase(0,1);
            }
            for(n=0; n<3; n++)
            {
               if(temp.count(atoms[n])==0)
               {
                  cout<<"Error 12amber "<<atoms[n]<<" "<<res<<endl;
                  exit(0);
               }
               sideatomindex[i].push_back(temp[atoms[n]]);
               sidechainpara[i].push_back(data[2-n]);
            }
            sideatom[i].push_back(atoms[3]);
            temp[atoms[3]]=j;
            j++;
         }
         else if(mark=="STOP")
         {
            if(i<residue_type_num)
            {
               //cout<<i<<"   "<<res<<endl;
               atom2index[i]=temp;
            }
            i++;
         }
         else
         {
            //cout<<mark<<" ";
         }
      }
      //The last word "STOP" will read twice! 
      if(i==residue_type_num+1) {}//cout<<"Top_rotamer.inp is read in..."<<endl;
      else {cout<<"Error in reading Top_rotamer.inp! "<<i<<endl; exit(0);}
      return 1;
}
Rotamer2Structure::Rotamer2Structure(const char* inp)
{
            
      vector<int> index; strs name;
      mainchainpara.assign(residue_type_num, xyz);
      sidechainpara.assign(residue_type_num, xyz);
      hydrogenspara.assign(residue_type_num, xyz);
      mainatomindex.assign(residue_type_num, index);
      sideatomindex.assign(residue_type_num, index);
      hydrogenindex.assign(residue_type_num, index);
      mainatom.assign(residue_type_num, name);
      sideatom.assign(residue_type_num, name);
      hydrogen.assign(residue_type_num, name);
      intercoor.assign(3, 0);
      cartcoor.assign(3, 0);
      amber_para=false;
      
      
      ifstream oi(inp, ios::in);
      if(!oi.good()) cout<<"Error in opening para file: "<<inp<<endl;
      string res, mark;
      strs atoms(4);
      str2int temp;
      vector<float> data(5);
      int i=0, j=0, n=0;
      float radii=0;
      while(oi.good())
      {
         oi>>mark;
         if(mark=="RESI")
         {
            oi>>res>>radii;
            res2index[res]=i;
            res2radii[res]=radii;
            index2res.push_back(res);
            //cout<<i<<"   "<<res<<endl;
            j=3;
            temp.clear();
            temp["-C"]=0; mainatom[i].push_back("-C");
            temp["N"]=1; mainatom[i].push_back("CA");
            temp["CA"]=2;  mainatom[i].push_back("N");
         }
#ifdef CB
         else if(mark=="IC" || mark=="CC")
#else
         else if(mark=="IC")
#endif
         {
            oi>>atoms[0]>>atoms[1]>>atoms[2]>>atoms[3]>>data[0]>>data[1]>>data[2]>>data[3]>>data[4];
            if(atoms[2][0]=='*')
            {
               atoms[2].erase(0,1);
            }
            for(n=0; n<3; n++)
            {
               if(temp.count(atoms[n])==0)
               {
                  cout<<"Error 11 "<<atoms[n]<<" "<<res<<endl;
                  exit(0);
               }
               mainatomindex[i].push_back(temp[atoms[n]]);
               mainchainpara[i].push_back(data[2+n]);
            }
            mainatom[i].push_back(atoms[3]);
            temp[atoms[3]]=j;
            j++;
         }
#ifdef CB
         else if(mark=="SC")
#else
         else if(mark=="SC" || mark=="CC")
#endif
         {
            oi>>atoms[0]>>atoms[1]>>atoms[2]>>atoms[3]
            >>data[0]>>data[1]>>data[2]>>data[3]>>data[4];
            if(atoms[2][0]=='*')
            {
               atoms[2].erase(0,1);
            }
            for(n=0; n<3; n++)
            {
               if(temp.count(atoms[n])==0)
               {
                  cout<<"Error 12 "<<atoms[n]<<" "<<res<<endl;
                  exit(0);
               }
               sideatomindex[i].push_back(temp[atoms[n]]);
               sidechainpara[i].push_back(data[2+n]);
            }
            sideatom[i].push_back(atoms[3]);
            temp[atoms[3]]=j;
            j++;
         }
         else if(mark=="HC")
         {
            oi>>atoms[0]>>atoms[1]>>atoms[2]>>atoms[3]
            >>data[0]>>data[1]>>data[2]>>data[3]>>data[4];
            if(atoms[2][0]=='*')
            {
               atoms[2].erase(0,1);
            }
            for(n=0; n<3; n++)
            {
               if(temp.count(atoms[n])==0)
               {
                  cout<<"Error 13 "<<atoms[n]<<" "<<res<<endl;
                  exit(0);
               }
               hydrogenindex[i].push_back(temp[atoms[n]]);
               hydrogenspara[i].push_back(data[2+n]);
            }
            hydrogen[i].push_back(atoms[3]);
            temp[atoms[3]]=j;
            j++;
         }
         else if(mark=="STOP")
         {
            if(i<residue_type_num)
            {
               //cout<<i<<"   "<<res<<endl;
               atom2index[i]=temp;
            }
            i++;
         }
         //else cout<<mark<<" ";
         
      }
      //The last word "STOP" will read twice! 
      if(i==residue_type_num+1) {}//cout<<"Top_rotamer.inp is read in..."<<endl;
      else {cout<<"Error in reading Top_rotamer.inp! "<<i<<endl; exit(0);}
      
}
//torsion contains Chi 
//For dee sidechain building!
float Rotamer2Structure::ConstructAASidechain(const string& res, str2vectord &NCACB, 
                       const vector <float> &torsion, str2vectord &cartesian)
{
#ifdef discal
      float dis=0, rec=0;
#endif
      if(res2index.count(res)==0)
      {
         cout<<"Rotamer2Structure::ConstructAASidechain  Error in residue type "<<res<<endl;
         exit(0);
      }
      int resi=res2index[res];
      int i=0, j=0, k=0;
      cycle.assign(sideatom[resi].size()+mainatom[resi].size(), xyz);
      cycle[1]=NCACB["N"];
      cycle[2]=NCACB["CA"];

#ifdef CB 
      cycle[8]=NCACB["CB"];
#else     
      cycle[4]=NCACB["C"];
#endif

      for(i=0; i<sideatom[resi].size(); i++)
      {
         j=i*3;
#ifdef CB 
         if( i<torsion.size() )
         {
            intercoor[2]=torsion[i];
         }
#else
         if( i-1<torsion.size() && i>0)
         {
            intercoor[2]=torsion[i-1];
         }
#endif
         else
         {
            intercoor[2]=sidechainpara[resi][j];
         }
         intercoor[1]=sidechainpara[resi][j+1];
         intercoor[0]=sidechainpara[resi][j+2];
         
         //Using Amber Default Parameters for different top relationships.
         if(amber_para && i==torsion.size())
         {
            switch(resi)
            {
               case 2:  intercoor[2]=torsion[i-1]+180; break; //ASN
               case 3:  intercoor[2]=torsion[i-1]+180; break; //ASP
               case 5:  intercoor[2]=torsion[i-1]+180; break; //GLN
               case 6:  intercoor[2]=torsion[i-1]+180; break; //GLU
               case 9:  intercoor[2]=torsion[i-2]-120; break; //ILE
               case 10: intercoor[2]=torsion[i-1]+120; break; //LEU
               case 16: intercoor[2]=torsion[i-1]+240; break; //THR
               case 19: intercoor[2]=torsion[i-1]+120; break; //VAL
               default: break;
            }
         }
         
         if(cycle[sideatomindex[resi][j]].size()!=3) cout<<"Error1 "<<sideatomindex[resi][j]<<"   "<<sideatom[resi][i]<<res<<endl;
         if(cycle[sideatomindex[resi][j+1]].size()!=3) cout<<"Error2 "<<sideatomindex[resi][j+1]<<"   "<<sideatom[resi][i]<<res<<endl;
         if(cycle[sideatomindex[resi][j+2]].size()!=3) cout<<"Error3 "<<sideatomindex[resi][j+2]<<"   "<<sideatom[resi][i]<<res<<endl;
         
         if(!internal2cartesian(             cycle[sideatomindex[resi][j]], 
            cycle[sideatomindex[resi][j+1]], cycle[sideatomindex[resi][j+2]], 
            intercoor, cartcoor))
         {
            cout<<i<<" "<<res<<" SEE "<<sideatomindex[resi][j]<<"  "<<sideatomindex[resi][j+1]<<"  "<<sideatomindex[resi][j+2]<<endl;
            cout<<"LOOK "<<sideatom[resi][i]<<" "<<mainatom[resi].size()+i<<endl;
            exit(0);
         }
         /*
         if(res=="ILE")
         {
            cout<<mainatom[resi].size()+i<<" ";
            ShowMyvector(intercoor);
            ShowMyvector(cartcoor);
            cout<<endl;
         }
         */
         cartesian[sideatom[resi][i]]=cartcoor;
         cycle[mainatom[resi].size()+i]=cartcoor;
#ifdef discal         
         dis=Points2Distance(cycle[8], cartcoor);
         if(dis>rec) {rec=dis;}
#endif
      }
#ifdef discal     
      return rec;
#else
      return  0;
#endif
}
float Rotamer2Structure::res2radius(string res)
{
      return res2radii[res];
}
int Rotamer2Structure::ConstructHN(const string& res, str2vectord &mC, 
                       str2vectord &NCA)
{ 
      //PRO is a different case for it has no HN but CD
      if(res=="PRO") return 0;
      
      string aa;
      if(res=="HIS") aa="HSD";
      else aa=res;
      if(res2index.count(aa)==0)
      {
         cout<<"Rotamer2Structure::ConstructHN Warning in residue "<<aa<<endl;
         return 0;
      }
      int resi=res2index[res];
      cycle[0]=mC["C"];
      cycle[1]=NCA["N"];
      cycle[2]=NCA["CA"];
            
      //HN calculation
      intercoor[2]=mainchainpara[resi][0];
      intercoor[1]=mainchainpara[resi][1];
      intercoor[0]=mainchainpara[resi][2];
      
      internal2cartesian(cycle[mainatomindex[resi][0]], cycle[mainatomindex[resi][1]], cycle[mainatomindex[resi][2]], intercoor, cartcoor);
      
      NCA[mainatom[resi][3]]=cartcoor;

      return 1;
}

int Rotamer2Structure::AtomIndexMaping(const string& res, str2vectord &coordata)
{
      if(res2index.count(res)==0)
      {
         cout<<"Rotamer2Structure::AtomIndexMaping Warning in residue type "<<res<<endl;
         return 0;
      }
      int resi=res2index[res]; str2vectord::iterator pos2;
      cycle.assign(sideatom[resi].size()+mainatom[resi].size()+hydrogen[resi].size(), xyz);
      for(pos2=coordata.begin(); pos2!=coordata.end(); ++pos2)
      {
         if(atom2index[resi].count(pos2->first)==0)
         {
            cout<<"Error in atom type in "<<res<<" "<<pos2->first<<endl;
            exit(0);
         }
         cycle[atom2index[resi][pos2->first]]=pos2->second;
      }
      return 1;
}

//Add H of side chain
int Rotamer2Structure::ConstructHs(const string& res, str2vectord &group) 
{
      if( AtomIndexMaping(res, group) == 0 ) return 0;
      int resi=res2index[res];
      int i=0, j=0;

      for(i=0; i<hydrogen[resi].size(); i++)
      {
         j=i*3;
         intercoor[2]=hydrogenspara[resi][j];
         intercoor[1]=hydrogenspara[resi][j+1];
         intercoor[0]=hydrogenspara[resi][j+2];

         if(cycle[hydrogenindex[resi][j]].size()!=3) 
            cout<<"Error4 "<<hydrogenindex[resi][j]<<"   "
            <<hydrogen[resi][i]<<res<<endl;
         if(cycle[hydrogenindex[resi][j+1]].size()!=3) 
            cout<<"Error5 "<<hydrogenindex[resi][j+1]<<"   "
            <<hydrogen[resi][i]<<res<<endl;
         if(cycle[hydrogenindex[resi][j+2]].size()!=3) 
            cout<<"Error6 "<<hydrogenindex[resi][j+2]<<"   "
            <<hydrogen[resi][i]<<res<<endl;
         
         if(!internal2cartesian(             cycle[hydrogenindex[resi][j]], 
            cycle[hydrogenindex[resi][j+1]], cycle[hydrogenindex[resi][j+2]], 
            intercoor, cartcoor))
         {
            cout<<i<<" "<<res<<" SEE "<<hydrogenindex[resi][j]<<"  "<<hydrogenindex[resi][j+1]<<"  "<<hydrogenindex[resi][j+2]<<endl;
            cout<<"LOOK "<<hydrogen[resi][i]<<" "<<mainatom[resi].size()+i<<endl;
            exit(0);
         }
         group[hydrogen[resi][i]]=cartcoor;
         cycle[atom2index[resi][hydrogen[resi][i]]]=cartcoor;
      }
      return i; //Number of hydrogens
}


int Rotamer2Structure::Atom2Output(str2strs &outseq)
{
      strs temp;
      for(int i=0; i<residue_type_num; ++i)
      {
         temp=mainatom[i];
         temp.insert(temp.end(), sideatom[i].begin(), sideatom[i].end());
         temp.insert(temp.end(), hydrogen[i].begin(), hydrogen[i].end());
         outseq[index2res[i]]=temp;
         temp.clear();
      }
      return 1;
}

int Rotamer2Structure::OutPdbFile(strs &peptide, array2sd &cartesian)
{
      int i, j, k, l;
      string swap;
      str2strs outseq;
      Atom2Output(outseq);
      str2vectord::iterator pos;
      for(i=0, k=0; i<cartesian.size(); ++i)
      {
         if(cartesian[i].size()>0 && outseq[peptide.at(i)].size()>0 )
         for(j=0; j<outseq[peptide.at(i)].size(); ++j)
         {
            //if(outseq[peptide.at(i)].at(j)[0]=='H') continue;
            if(cartesian.at(i).count(outseq[peptide.at(i)].at(j))!=1)
            continue;
            if(peptide.at(i)=="PRO" && j==3) continue;
            
            if(peptide.at(i)=="HSD")
            PDB_Format_OUT(k+1, outseq[peptide.at(i)].at(j), "HIS", i+1, (cartesian.at(i)[outseq[peptide.at(i)].at(j)]) );
            else if(peptide.at(i)=="LEU" && outseq[peptide.at(i)].at(j)=="CD1")
            {
               if(cartesian.at(i).count("CD2")!=1) continue;
               PDB_Format_OUT(k+1, outseq[peptide.at(i)].at(j), peptide.at(i), i+1, (cartesian.at(i)["CD2"]) );
            }
            else if(peptide.at(i)=="LEU" && outseq[peptide.at(i)].at(j)=="CD2")
            {
               if(cartesian.at(i).count("CD1")!=1) continue;
               PDB_Format_OUT(k+1, outseq[peptide.at(i)].at(j), peptide.at(i), i+1, (cartesian.at(i)["CD1"]) );
            }
            else if(peptide.at(i)=="ILE" && outseq[peptide.at(i)].at(j)=="CD")
            {
               PDB_Format_OUT(k+1, "CD1", peptide.at(i), i+1, (cartesian.at(i)[outseq[peptide.at(i)].at(j)]) );
            }
            else if(outseq[peptide.at(i)].at(j)[0]=='H'
            && outseq[peptide.at(i)].at(j).size()==4)
            {
               swap=outseq[peptide.at(i)].at(j);
               swap[0]=outseq[peptide.at(i)].at(j)[3];
               for(l=0; l<3; l++)
                  swap[l+1]=outseq[peptide.at(i)].at(j)[l];
               PDB_Format_OUT(k+1, swap, peptide.at(i), i+1, (cartesian.at(i)[outseq[peptide.at(i)].at(j)]) );
            }
            else
            PDB_Format_OUT(k+1, outseq[peptide.at(i)].at(j), peptide.at(i), i+1, (cartesian.at(i)[outseq[peptide.at(i)].at(j)]) );
            ++k;
         }
         else
         {
            for(pos=cartesian.at(i).begin(); pos!=cartesian.at(i).end(); ++pos)
            {
               PDB_Format_OUT(k+1, pos->first, peptide.at(i), i+1, pos->second);
               ++k;
            }
         }
      }
      return 1;
}

int Rotamer2Structure::OutPdbFile(strs &peptide, array2sd &cartesian, ofstream &pdbout)
{
      int i, j, k, l;
      string swap;
      str2strs outseq;
      Atom2Output(outseq);
      str2vectord::iterator pos;
      for(i=0, k=0; i<cartesian.size(); ++i)
      {
         if(cartesian[i].size()>0 && outseq[peptide.at(i)].size()>0 )
         for(j=0; j<outseq[peptide.at(i)].size(); ++j)
         {
            //if(outseq[peptide.at(i)].at(j)[0]=='H') continue;
            if(cartesian.at(i).count(outseq[peptide.at(i)].at(j))!=1)
            continue;
            
            if(peptide.at(i)=="PRO" && j==3) continue; 
            
            if(peptide.at(i)=="HSD")
            PDB_Format_OUT(k+1, outseq[peptide.at(i)].at(j), "HIS", i+1, (cartesian.at(i)[outseq[peptide.at(i)].at(j)]), pdbout);
            else if(peptide.at(i)=="LEU" && outseq[peptide.at(i)].at(j)=="CD1")
            {
               if(cartesian.at(i).count("CD2")!=1) continue;
               PDB_Format_OUT(k+1, outseq[peptide.at(i)].at(j), peptide.at(i), i+1, (cartesian.at(i)["CD2"]) , pdbout);
            }
            else if(peptide.at(i)=="LEU" && outseq[peptide.at(i)].at(j)=="CD2")
            {
               if(cartesian.at(i).count("CD1")!=1) continue;
               PDB_Format_OUT(k+1, outseq[peptide.at(i)].at(j), peptide.at(i), i+1, (cartesian.at(i)["CD1"]) , pdbout);
            }
            else if(peptide.at(i)=="ILE" && outseq[peptide.at(i)].at(j)=="CD")
            {
               PDB_Format_OUT(k+1, "CD1", peptide.at(i), i+1, (cartesian.at(i)[outseq[peptide.at(i)].at(j)]) , pdbout);
            }
            else if(outseq[peptide.at(i)].at(j)[0]=='H'
            && outseq[peptide.at(i)].at(j).size()==4)
            {
               swap=outseq[peptide.at(i)].at(j);
               swap[0]=outseq[peptide.at(i)].at(j)[3];
               for(l=0; l<3; l++)
                  swap[l+1]=outseq[peptide.at(i)].at(j)[l];
               PDB_Format_OUT(k+1, swap, peptide.at(i), i+1, (cartesian.at(i)[outseq[peptide.at(i)].at(j)]) , pdbout);
            }
            else if(outseq[peptide.at(i)].at(j)=="HN")
            {
               PDB_Format_OUT(k+1, "H", peptide.at(i), i+1, (cartesian.at(i)[outseq[peptide.at(i)].at(j)]) , pdbout);
            }
            else
            PDB_Format_OUT(k+1, outseq[peptide.at(i)].at(j), peptide.at(i), i+1, (cartesian.at(i)[outseq[peptide.at(i)].at(j)]), pdbout);
            ++k;
         }
         else
         {
            for(pos=cartesian.at(i).begin(); pos!=cartesian.at(i).end(); ++pos)
            {
               PDB_Format_OUT(k+1, pos->first, peptide.at(i), i+1, pos->second, pdbout);
               ++k;
            }
         }
      }
      return 1;
}

int Rotamer2Structure::OutPdbFile(strs &peptide, array2sd &cartesian, vector<char> &chn, ofstream &pdbout)
{
      int i, j, k, l;
      string swap;
      str2strs outseq;
      Atom2Output(outseq);
      str2vectord::iterator pos;
      for(i=0, k=0; i<cartesian.size(); ++i)
      {
         if(cartesian[i].size()>0 && outseq[peptide.at(i)].size()>0 )
         for(j=0; j<outseq[peptide.at(i)].size(); ++j)
         {
            //if(outseq[peptide.at(i)].at(j)[0]=='H') continue;
            if(cartesian.at(i).count(outseq[peptide.at(i)].at(j))!=1)
            continue;
            
            if(peptide.at(i)=="PRO" && j==3) continue; 
            
            if(peptide.at(i)=="HSD")
            PDB_Format_OUT(k+1, outseq[peptide.at(i)].at(j), "HIS", chn[i], i+1, (cartesian.at(i)[outseq[peptide.at(i)].at(j)]), pdbout);
            else if(peptide.at(i)=="LEU" && outseq[peptide.at(i)].at(j)=="CD1")
            {
               if(cartesian.at(i).count("CD2")!=1) continue;
               PDB_Format_OUT(k+1, outseq[peptide.at(i)].at(j), peptide.at(i), chn[i], i+1, (cartesian.at(i)["CD2"]) , pdbout);
            }
            else if(peptide.at(i)=="LEU" && outseq[peptide.at(i)].at(j)=="CD2")
            {
               if(cartesian.at(i).count("CD1")!=1) continue;
               PDB_Format_OUT(k+1, outseq[peptide.at(i)].at(j), peptide.at(i), chn[i], i+1, (cartesian.at(i)["CD1"]) , pdbout);
            }
            else if(peptide.at(i)=="ILE" && outseq[peptide.at(i)].at(j)=="CD")
            {
               PDB_Format_OUT(k+1, "CD1", peptide.at(i),chn[i], i+1, (cartesian.at(i)[outseq[peptide.at(i)].at(j)]) , pdbout);
            }
            else if(outseq[peptide.at(i)].at(j)[0]=='H'
            && outseq[peptide.at(i)].at(j).size()==4)
            {
               swap=outseq[peptide.at(i)].at(j);
               swap[0]=outseq[peptide.at(i)].at(j)[3];
               for(l=0; l<3; l++)
                  swap[l+1]=outseq[peptide.at(i)].at(j)[l];
               PDB_Format_OUT(k+1, swap, peptide.at(i), chn[i], i+1, (cartesian.at(i)[outseq[peptide.at(i)].at(j)]) , pdbout);
            }
            else if(outseq[peptide.at(i)].at(j)=="HN")
            {
               PDB_Format_OUT(k+1, "H", peptide.at(i), chn[i], i+1, (cartesian.at(i)[outseq[peptide.at(i)].at(j)]) , pdbout);
            }
            else
            PDB_Format_OUT(k+1, outseq[peptide.at(i)].at(j), peptide.at(i), chn[i], i+1, (cartesian.at(i)[outseq[peptide.at(i)].at(j)]), pdbout);
            ++k;
         }
         else
         {
            for(pos=cartesian.at(i).begin(); pos!=cartesian.at(i).end(); ++pos)
            {
               PDB_Format_OUT(k+1, pos->first, peptide.at(i), chn[i], i+1, pos->second, pdbout);
               ++k;
            }
         }
      }
      return 1;
}

int Rotamer2Structure::OutPdbFileALL(strs &peptide, array2sd &cartesian, ofstream &pdbout)
{
      int i, j, k, l;
      string swap;
      str2strs outseq;
      Atom2Output(outseq);
      str2vectord::iterator pos;
      for(i=0, k=0; i<cartesian.size(); ++i)
      {
         for(pos=cartesian.at(i).begin(); pos!=cartesian.at(i).end(); ++pos)
         {
            PDB_Format_OUT(k+1, pos->first, peptide.at(i), i+1, pos->second, pdbout);
            ++k;
         }
      }
      return 1;
}

int Rotamer2Structure::OutPdbFileALL(strs &peptide, array2sd &cartesian, vector<char> &chn, ofstream &pdbout)
{
      int i, j, k, l;
      string swap;
      str2strs outseq;
      Atom2Output(outseq);
      str2vectord::iterator pos;
      for(i=0, k=0; i<cartesian.size(); ++i)
      {
         for(pos=cartesian.at(i).begin(); pos!=cartesian.at(i).end(); ++pos)
         {
            PDB_Format_OUT(k+1, pos->first, peptide.at(i), chn[i], i+1, pos->second, pdbout);
            ++k;
         }
      }
      return 1;
}
#endif




