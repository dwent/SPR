/*************************************************************************
 in  GetData(char *input) modified to detect number or char
 if(strcmp(temp, " ")!=0        => 
 if(strcmp(temp, " ")!=0 && (int)temp[0]>47 && (int)temp[0]<58)
            strcat(rac[i].atom, temp);
            else if(strcmp(temp, " ")!=0 && (int)temp[0]>64 && (int)temp[0]<91)
            {
               initialization(rac[i].atom, 0, 7);
               copy(rac[i].atom, line, 12, 16);
            }
 Date 2007.5.10
 
 put racfree() in GetData(char *input); GetDecoyData(char *input, char* file2);
 add RemoveHydrogen();
 Date 2007.5.17. 
 
 in VerifyBackboneData() add else { return 0; }
 Date 2007.6.19.
 CB belong to main chain atom
 Date 2007.11.2.
 substitute residue in struct rac for protein design
 Date 2008.3.3.
 Modefied algorithm of phi psi in
 int ProteinStructureAnalysis::PhiPsi(...)
 for rubust to the sequence of main chain atoms.
 
 增加GetData2，用于对Charmm命名方式的PDB文件读取
 Date: 2008.8.1.
 增加PhiPsiPro2, 用于对多条链的蛋白phipsiomg计算
 但测试发现原先的PhiPsiPro也能算，而且不受残基序号重复的影响，也能算好。
 ExtractPDB_Gapless用于提取无间隔的序列与结构数据
 Date: 2008.8.19
 取消了CalRotamer相关部分
 2008.10.27.
 GetData2少了x1=NULL; atom_num=i; 补上
 Date2008.11.3. 
 简化这个类，不需要Geometry, CalRotamer等（本来就用的少）
 Date: 2008.11.9.
 增加：
 int ExtractPDB(strs &peptide, array2sd &cartesian, myvecint &seq_index);
 int ExtractChain(vector<char> chn);//给出每个残基的链号
 Date: 2008.12.17.
 发现BUG： 在GetData2中没有rac[i].chain=line[21];加上
 RemoveSideChain()和RemoveHydrogen()中关于chain的修改，原先是赋值‘ ’
 Date: 2008.12.23.
**************************************************************************/

#ifndef ProteinToolsHead
#define ProteinToolsHead


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
using namespace std;


const int Max_Atom_Number_BS=100000;

#ifdef CB
const bool CBinclude=true;
#else
const bool CBinclude=false;
#endif
/******************************************************************************
class ProteinStructureAnalysis 
Check PDB file, Obtain Sequence and Coordinates
Calculate Phi, Psi, Omega dihedral angles
******************************************************************************/
class ProteinStructureAnalysis
{
      public: 
      ProteinStructureAnalysis();
      ~ProteinStructureAnalysis();

      void initialization(char *str, int l1, int l2);
      void copy(char *str1, char *str2, int x, int y);
      void XYZ(int i, string &atom, vector<float> &xyz);
      void XYZ(int i, string &atom, string &resi, vector<float> &xyz);

      void racfree();
      char Chain(int i);
      
      int Atom_Num();
      int Seq(int i);
      int GetData(const char *input);
      int GetData2(const char *input);
      int GetDecoyData(char *input, char* file2);
      int TransferData(const strs *residue_name, array2sd *cartesian);
      int Coordinates(array2sd *cartesian);
      int ExtractPDB(strs *peptide, array2sd *cartesian);
      int ExtractPDB_Gapless(strs *peptide, array2sd *cartesian);
      int ExtractPDB(strs &peptide, array2sd &cartesian, myvecint &seq_index);
      int ExtractChain(vector<char> &chn);//给出每个残基的链号
      int RemoveHydrogen(array2sd &cartesian);
      int RemoveHydrogen();
      int RemoveSideChain();
      
      int VerifyPdbData();
      int VerifyPdbData(vector<int> &States);
      int VerifyBackboneData();
      int ReadProteinBackbone(array2sd *BackboneAtom); //不识别多链情况，做接在一起处理
      int ReadProteinSequence(strs *rps); //不识别多链情况，做接在一起处理
      int ReadProteinSequence(strs *rps, vector<int>* index); //不识别多链情况，做接在一起处理
      int SubstituteResidue(const vector<int> &index, const strs &substitution);
      
      protected:
      
      ResiAtomCoor_BS *rac;
      int atom_num;
};



/*****************************************************************************************
   int ProteinStructureAnalysis::ExtractPDB(strs *peptide, array2sd *cartesian)
   ADD:   peptide->clear(); cartesian->clear();
   Date: 2008.4.1.
   ADD:   GetData(): LEU CD1->CD2 and CD2->CD1 
   This is because Naming system difference. Refer: Psfgen
   Date: 2008.4.24.
   Bug in ReadProteinSequence(strs *rps, vector<int>* index)
   ExtractPDB 增加了对序号不连续的序列的多肽的处理，当出现不连续时候，peptide, 
   cartesian增加一个空的元素。
   Date: 2008.5.10.
   增加GetData2，用于对Charmm命名方式的PDB文件读取
   Date: 2008.8.1.
   GetData增加OT1,OT2处理。
   Date2008.9.8.
******************************************************************************************/

ProteinStructureAnalysis::ProteinStructureAnalysis()
{
      int i=0;
      rac=new ResiAtomCoor_BS[Max_Atom_Number_BS];
      for(i=0; i<Max_Atom_Number_BS; i++)
      {
         rac[i].sequence=0;
         rac[i].chain=' ';
         initialization(rac[i].resi, 0, 6);
         initialization(rac[i].atom, 0, 7);
         rac[i].xyz[0]=0; rac[i].xyz[1]=0; rac[i].xyz[2]=0;
      }
      atom_num=0;
}

void ProteinStructureAnalysis::racfree()
{
      int i=0; atom_num=0;
      for(i=0; i<Max_Atom_Number_BS; i++)
      {
         rac[i].sequence=0;
         rac[i].chain=' ';
         initialization(rac[i].resi, 0, 6);
         initialization(rac[i].atom, 0, 7);
         rac[i].xyz[0]=0; rac[i].xyz[1]=0; rac[i].xyz[2]=0;
      }
}

ProteinStructureAnalysis::~ProteinStructureAnalysis()
{
      delete[] rac;
}

void ProteinStructureAnalysis::initialization(char *str, int l1, int l2)
{
      for(int i=l1; i<l2; i++)
      {
         str[i]='\0';
      }
}

void ProteinStructureAnalysis::copy(char *str1, char *str2, int x, int y)
{
      int i=0;
      int j=0;
      for(i=x; i<y; i++)
      {
         if(str2[i]!=' ')
         {
            str1[j]=str2[i];
            //cout<<str1[j]<<"|";
            j++;
         }
      }
      for(i=j; i<y-x; i++)
      {
         str1[i]='\0';
      }
}


int ProteinStructureAnalysis::GetData(const char *input)
{
      int i=0;
      char line[100], temp[8];
      char *x1;
      
      ifstream oi (input, ios::in);
      if(!oi)
      {
         cout<<"Error in file reading "<<input<<endl;
         return 0;
      }
      racfree();
      oi.seekg(0,ios::beg);
      while(oi.good())
      {
         oi.getline(line, 100, '\n');   
         x1=strstr(line, "ATOM");
         if(x1 != NULL && line[0]=='A')
         {
            if(strstr(line, "OXT")!=0) continue;
            if(strstr(line, "OT2")!=0) continue;
            
            copy(rac[i].resi, line, 17, 20);
            copy(rac[i].atom, line, 13, 16);
            
            initialization(temp, 0, 8);
            copy(temp, line, 12, 13);
            if(strcmp(temp, " ")!=0 && (int)temp[0]>47 && (int)temp[0]<58)
            strcat(rac[i].atom, temp);
            else if(strcmp(temp, " ")!=0 && (int)temp[0]>64 && (int)temp[0]<91)
            {
               initialization(rac[i].atom, 0, 7);
               copy(rac[i].atom, line, 12, 16);
            }
            
            if(strcmp(rac[i].atom, "H")==0) 
            {
               //cout<<"H:"<<rac[i].atom<<" "<<rac[i].resi<<"\n";
               strcat(rac[i].atom, "N");
            }
            if(strcmp(rac[i].atom, "HG")==0 && strcmp(rac[i].resi, "SER")==0)
            {
               //cout<<"HG:"<<rac[i].atom<<" "<<rac[i].resi<<"\n";
               strcat(rac[i].atom, "1");
            }
            if(strcmp(rac[i].atom, "CD1")==0 && strcmp(rac[i].resi, "ILE")==0)
            {
               rac[i].atom[2]= '\0';
            }
            if(strcmp(rac[i].atom, "CD1")==0 && strcmp(rac[i].resi, "LEU")==0)
            {
               rac[i].atom[2]= '2';
            }
            else if(strcmp(rac[i].atom, "CD2")==0 && strcmp(rac[i].resi, "LEU")==0)
            {
               rac[i].atom[2]= '1';
            }
            if(strcmp(rac[i].resi, "HIS")==0)
            {
               initialization(rac[i].resi, 0, 6);
               strcat(rac[i].resi, "HSD");
            }
            if(strcmp(rac[i].resi, "ILE")==0)
            {
               if(strcmp(rac[i].atom, "HD11")==0)
               {
                  rac[i].atom[3]='\0';
               }
               else if(strcmp(rac[i].atom, "HD12")==0)
               {
                  rac[i].atom[3]='\0';
                  rac[i].atom[2]='2';
               }
               else if(strcmp(rac[i].atom, "HD13")==0 )
               {
                  rac[i].atom[3]='\0';
                  rac[i].atom[2]='3';
               }
            }
            if(strcmp(rac[i].atom, "OT1")==0) 
            {
               rac[i].atom[1]='\0';
               rac[i].atom[2]='\0';
               rac[i].atom[3]='\0';
            }
            
            
            rac[i].chain=line[21];
            
            initialization(temp, 0, 8);
            copy(temp, line, 22, 26);
            rac[i].sequence=atoi(temp);
            
            initialization(temp, 0, 8);
            copy(temp, line, 30, 38);
            rac[i].xyz[0]=atof(temp);
            
            initialization(temp, 0, 8);
            copy(temp, line, 38, 46);
            rac[i].xyz[1]=atof(temp);
            
            initialization(temp, 0, 8);
            copy(temp, line, 46, 54);
            rac[i].xyz[2]=atof(temp);
            
            i++;
            if( i>Max_Atom_Number_BS )
            {
               cout<<"Error in GetData(). Overflow from the Max_Atom_Number_BS limit!\n";
               break;
            }
         }
      }
      /*Verify the data.
      for(int j=0; j<i; j++)
      {
         cout<<rac[j].sequence<<" "<<rac[j].resi<<" "
         <<rac[j].atom<<" "<<rac[j].xyz[0]<<" "<<rac[j].xyz[1]<<" "<<rac[j].xyz[2]<<"\n";
      }
      //*/
      x1=NULL;
      atom_num=i;
      return i;
}

int ProteinStructureAnalysis::GetData2(const char *input)
{
      int i=0;
      char line[100], temp[8];
      char *x1;
      
      ifstream oi (input, ios::in);
      if(!oi)
      {
         cout<<"Error in file reading "<<input<<endl;
         exit(0);
      }
      oi.seekg(0,ios::beg);
      while(oi.good())
      {
         oi.getline(line, 100, '\n');   
         x1=strstr(line, "ATOM");
         if(x1 != NULL && line[0]=='A')
         {
            copy(rac[i].resi, line, 17, 20);
            copy(rac[i].atom, line, 12, 16);
            
            initialization(temp, 0, 8);
            copy(temp, line, 22, 26);
            rac[i].sequence=atoi(temp);
            
            initialization(temp, 0, 8);
            copy(temp, line, 30, 38);
            rac[i].xyz[0]=atof(temp);
            
            initialization(temp, 0, 8);
            copy(temp, line, 38, 46);
            rac[i].xyz[1]=atof(temp);
            
            initialization(temp, 0, 8);
            copy(temp, line, 46, 54);
            rac[i].xyz[2]=atof(temp);
            
            rac[i].chain=line[21];
            
            i++;
            if( i>Max_Atom_Number_BS )
            {
               cout<<"Error in GetData2(). Overflow from the Max_Atom_Number limit!\n";
            }
         }
      }
      /*Verify the data.
      for(int j=0; j<i; j++)
      {
         cout<<rac[j].sequence<<" "<<rac[j].resi<<" "
         <<rac[j].atom<<" "<<rac[j].xyz[0]<<" "<<rac[j].xyz[1]<<" "<<rac[j].xyz[2]<<"\n";
      }
      */
      x1=NULL;
      atom_num=i;
      return i;
}

int ProteinStructureAnalysis::ExtractPDB(strs *peptide, array2sd *cartesian)
{
      vector<float> coors;
      str2vectord tmp_acid;
      int i=0, j=0;
      peptide->clear(); cartesian->clear();
      
      while(strlen(rac[i].resi)>0)
      {
         if(i==0)
         {
            string residue(rac[i].resi);
            peptide->push_back(residue);
            j=rac[i].sequence;
         }
         if(j+1==rac[i].sequence)
         {
            string residue(rac[i].resi);
            peptide->push_back(residue);
            j=rac[i].sequence;
            
            cartesian->push_back(tmp_acid);
            tmp_acid.clear();
         }
         else if(j!=rac[i].sequence)
         {
            string residuenull;
            peptide->push_back(residuenull);
            cartesian->push_back(tmp_acid);
            tmp_acid.clear();
            
            string residue(rac[i].resi);
            peptide->push_back(residue);
            j=rac[i].sequence;
            cartesian->push_back(tmp_acid);
         }
         string atom(rac[i].atom);
         coors.push_back(rac[i].xyz[0]);
         coors.push_back(rac[i].xyz[1]);
         coors.push_back(rac[i].xyz[2]);
         
         tmp_acid[atom]=coors;
         coors.clear();
         
         ++i;
      }
      if(strlen(rac[0].resi)>0)
         cartesian->push_back(tmp_acid);
      
      return i;
}


int ProteinStructureAnalysis::ExtractPDB_Gapless(strs *peptide, array2sd *cartesian)
{
      vector<float> coors;
      str2vectord tmp_acid;
      int i=0, j=0;
      peptide->clear(); cartesian->clear();
      
      while(strlen(rac[i].resi)>0)
      {
         if(i==0)
         {
            string residue(rac[i].resi);
            peptide->push_back(residue);
            j=rac[i].sequence;
         }
         if(j+1==rac[i].sequence)
         {
            string residue(rac[i].resi);
            peptide->push_back(residue);
            j=rac[i].sequence;
            
            cartesian->push_back(tmp_acid);
            tmp_acid.clear();
         }
         else if(j!=rac[i].sequence)
         {
            cout<<"GAP "<<j<<" "<<rac[i].sequence<<endl;
            string residue(rac[i].resi);
            peptide->push_back(residue);
            j=rac[i].sequence;
            cartesian->push_back(tmp_acid);
            tmp_acid.clear();
         }
         string atom(rac[i].atom);
         coors.push_back(rac[i].xyz[0]);
         coors.push_back(rac[i].xyz[1]);
         coors.push_back(rac[i].xyz[2]);
         
         tmp_acid[atom]=coors;
         coors.clear();
         
         ++i;
      }
      if(strlen(rac[0].resi)>0)
         cartesian->push_back(tmp_acid);
      
      return i;
}

int ProteinStructureAnalysis::ExtractPDB(strs &peptide, array2sd &cartesian, myvecint &seq_index)
{
      vector<float> coors;
      str2vectord tmp_acid;
      int i=0, j=0;
      peptide.clear(); cartesian.clear();
      
      while(strlen(rac[i].resi)>0)
      {
         if(i==0)
         {
            string residue(rac[i].resi);
            peptide.push_back(residue);
            seq_index.push_back(rac[i].sequence);
            j=rac[i].sequence;
         }
         if(j+1==rac[i].sequence)
         {
            string residue(rac[i].resi);
            peptide.push_back(residue);
            seq_index.push_back(rac[i].sequence);
            j=rac[i].sequence;
            
            cartesian.push_back(tmp_acid);
            tmp_acid.clear();
         }
         else if(j!=rac[i].sequence)
         {
            string residuenull;
            peptide.push_back(residuenull);
            seq_index.push_back(0);
            cartesian.push_back(tmp_acid);
            tmp_acid.clear();
            
            string residue(rac[i].resi);
            peptide.push_back(residue);
            seq_index.push_back(rac[i].sequence);
            j=rac[i].sequence;
            cartesian.push_back(tmp_acid);
         }
         string atom(rac[i].atom);
         coors.push_back(rac[i].xyz[0]);
         coors.push_back(rac[i].xyz[1]);
         coors.push_back(rac[i].xyz[2]);
         
         tmp_acid[atom]=coors;
         coors.clear();
         
         ++i;
      }
      if(strlen(rac[0].resi)>0)
         cartesian.push_back(tmp_acid);
      
      return i;
}
int ProteinStructureAnalysis::ExtractChain(vector<char> &chn)
{
      int i=0, j=0;
      chn.clear();
      
      while(strlen(rac[i].resi)>0)
      {
         if(i==0)
         {
            j=rac[i].sequence;
            chn.push_back(Chain(i));
         }
         if(j+1==rac[i].sequence)
         {
            j=rac[i].sequence;
            chn.push_back(Chain(i));
         }
         else if(j!=rac[i].sequence)
         {
            j=rac[i].sequence;
            chn.push_back(Chain(i));
            chn.push_back(Chain(i));
         }
         ++i;
      }
      return i; 
}


int ProteinStructureAnalysis::VerifyPdbData()
{
      int i=0, j=0, temp=0, n=0;
      int backbone[3]={0, 0, 0};
      
      typedef map<string, int> str2int;
      str2int resnum;
      resnum["ALA"]=5;
      resnum["ARG"]=11;
      resnum["ASN"]=8;
      resnum["ASP"]=8;
      resnum["CYS"]=6;
      resnum["GLN"]=9;
      resnum["GLU"]=9;
      resnum["GLY"]=4;
      resnum["HSD"]=10;
      resnum["HSE"]=10;
      resnum["HSP"]=10;
      resnum["ILE"]=8;
      resnum["LEU"]=8;
      resnum["LYS"]=9;
      resnum["MET"]=8;
      resnum["PHE"]=11;
      resnum["PRO"]=7;
      resnum["SER"]=6;
      resnum["THR"]=7;
      resnum["TRP"]=14;
      resnum["TYR"]=12;
      resnum["VAL"]=7;
      resnum["HIS"]=10;
      
      while(strlen(rac[i].resi)>0)
      {
         //cout<<i<<" "<<rac[i].resi<<" "<<rac[i].atom<<endl;
         if( rac[i].sequence==temp || i==0 )
         {
            if( strcmp(rac[i].atom, "N")==0) backbone[0]=1;
            if( strcmp(rac[i].atom, "C")==0) backbone[1]=1;
            if( strcmp(rac[i].atom, "CA")==0) backbone[2]=1;
            temp=rac[i].sequence;
            
            if(rac[i].atom[0]!='H')
               ++n;
            //cout<<"E1 "<<backbone[0]<<"*"<<backbone[1]<<"*"<<backbone[2]<<endl;
         }
         else if(rac[i].sequence==temp+1)
         {
            if(backbone[0]!=1 || backbone[1]!=1 || backbone[2]!=1)
            {
               cout<<"Error Format! (Lack of Backbone Atom)"<<endl;
               cout<<rac[i].sequence<<" "<<rac[i].resi<<" "<<rac[i].atom<<endl;
               return 0;
            }
            for(j=0; j<3; ++j)
               backbone[j]=0;
            
            //cout<<rac[i-1].resi<<" "<<"N= "<<n<<endl;
            str2int::iterator pos;
            pos=resnum.find(rac[i-1].resi);
            if(pos->second != n)
            {
               cout<<"Error! Lack of Atom at this residue "<<rac[i-1].resi<<endl;
               cout<<"Standard No. "<<resnum[rac[i-1].resi]<<" But it is "<<n<<endl;
               return 0;
            }
            
            n=0;
            
            if( strcmp(rac[i].atom, "N")==0) backbone[0]=1;
            if( strcmp(rac[i].atom, "C")==0) backbone[1]=1;
            if( strcmp(rac[i].atom, "CA")==0) backbone[2]=1;
            temp=rac[i].sequence;
            if(rac[i].atom[0]!='H')
               ++n;
         }
         else
         {
            cout<<"VerifyPdbData() "<<i<<" "<<rac[i].sequence<<"<->"<<temp<<" "<<rac[i].resi<<endl;
            cout<<"Error Format! (Noncontinous)"<<endl;
            return 0;
         }
         ++i;
      }
      
      str2int::iterator pos;
      pos=resnum.find(rac[i-1].resi);
      if(pos->second != n)
      {
         cout<<"Error! Lack of Atom at this residue "<<rac[i-1].resi<<endl;
         cout<<"Standard No. "<<resnum[rac[i-1].resi]<<" But it is "<<n<<endl;
         return 0;
      }
      
      //cout<<"OK! Confirmed PDB Format Protein Structure!"<<endl;
      return 1;
}
int ProteinStructureAnalysis::VerifyPdbData(vector<int> &States)
{
      int i=0, j=0, temp=0, n=0;
      int backbone[3]={0, 0, 0};
      
      typedef map<string, int> str2int;
      str2int resnum;
      resnum["ALA"]=5;
      resnum["ARG"]=11;
      resnum["ASN"]=8;
      resnum["ASP"]=8;
      resnum["CYS"]=6;
      resnum["GLN"]=9;
      resnum["GLU"]=9;
      resnum["GLY"]=4;
      resnum["HSD"]=10;
      resnum["HSE"]=10;
      resnum["HSP"]=10;
      resnum["ILE"]=8;
      resnum["LEU"]=8;
      resnum["LYS"]=9;
      resnum["MET"]=8;
      resnum["PHE"]=11;
      resnum["PRO"]=7;
      resnum["SER"]=6;
      resnum["THR"]=7;
      resnum["TRP"]=14;
      resnum["TYR"]=12;
      resnum["VAL"]=7;
      resnum["HIS"]=10;
      
      while(strlen(rac[i].resi)>0)
      {
         //cout<<i<<" "<<rac[i].resi<<" "<<rac[i].atom<<endl;
         if( rac[i].sequence==temp || i==0 )
         {
            if( strcmp(rac[i].atom, "N")==0) backbone[0]=1;
            if( strcmp(rac[i].atom, "C")==0) backbone[1]=1;
            if( strcmp(rac[i].atom, "CA")==0) backbone[2]=1;
            temp=rac[i].sequence;
            
            if(rac[i].atom[0]!='H' && strcmp(rac[i].atom, "OXT")!=0)
               ++n;
            //cout<<"E1 "<<backbone[0]<<"*"<<backbone[1]<<"*"<<backbone[2]<<endl;
         }
         else if(rac[i].sequence==temp+1)
         {
            if(backbone[0]!=1 || backbone[1]!=1 || backbone[2]!=1)
            {
               cout<<"Error Format! (Lack of Backbone Atom)"<<endl;
               cout<<rac[i].sequence<<" "<<rac[i].resi<<" "<<rac[i].atom<<endl;
               return 0;
            }
            for(j=0; j<3; ++j)
               backbone[j]=0;
            
            //cout<<rac[i-1].resi<<" "<<"N= "<<n<<endl;
            str2int::iterator pos;
            pos=resnum.find(rac[i-1].resi);
            if(pos->second != n)
            {
               cout<<"Error! Lack of Atom at this residue "<<rac[i-1].sequence<<rac[i-1].resi<<endl;
               cout<<"Standard No. "<<resnum[rac[i-1].resi]<<" But it is "<<n<<endl;
               States.push_back(0);
            }
            else
            {
               States.push_back(1);
            }
            n=0;
            
            if( strcmp(rac[i].atom, "N")==0) backbone[0]=1;
            if( strcmp(rac[i].atom, "C")==0) backbone[1]=1;
            if( strcmp(rac[i].atom, "CA")==0) backbone[2]=1;
            temp=rac[i].sequence;
            if(rac[i].atom[0]!='H' && strcmp(rac[i].atom, "OXT")!=0)
               ++n;
            
         }
         else
         {
            cout<<"VerifyPdbData() "<<i<<" "<<rac[i].sequence<<"<->"<<temp<<" "<<rac[i].resi<<endl;
            cout<<"Error Format! (Noncontinous)"<<endl;
            return 0;
         }
         ++i;
      }
      
      str2int::iterator pos;
      pos=resnum.find(rac[i-1].resi);
      if(pos->second != n)
      {
         cout<<"Error! Lack of Atom at this residue "<<rac[i-1].resi<<endl;
         cout<<"Standard No. "<<resnum[rac[i-1].resi]<<" But it is "<<n<<endl;
         States.push_back(0);
      }
      else
      {
         States.push_back(1);
      }
      //cout<<"OK! Confirmed PDB Format Protein Structure!"<<endl;
      return 1;
}


int ProteinStructureAnalysis::VerifyBackboneData()
{
      int i=0, j=0, temp=0, n=0;
      int backbone[4]={0, 0, 0, 0};
      
      while(strlen(rac[i].resi)>0)
      {
         //cout<<i<<" "<<rac[i].resi<<" "<<rac[i].atom<<endl;
         if( rac[i].sequence==temp || i==0 )
         {
            if( strcmp(rac[i].atom, "N")==0) backbone[0]=1;
            if( strcmp(rac[i].atom, "C")==0) backbone[1]=1;
            if( strcmp(rac[i].atom, "CA")==0) backbone[2]=1;
            if( strcmp(rac[i].atom, "O")==0) backbone[3]=1;
            temp=rac[i].sequence;
            //cout<<"E1 "<<backbone[0]<<"*"<<backbone[1]<<"*"<<backbone[2]<<endl;
         }
         else if(rac[i].sequence==temp+1)
         {
            if(backbone[0]!=1 || backbone[1]!=1 || backbone[2]!=1 || backbone[3]!=1)
            {
               cout<<"Error Format! (Lack of Backbone Atom)"<<endl;
               cout<<rac[i].sequence<<" "<<rac[i].resi<<" "<<rac[i].atom<<endl;
               cout<<backbone[0]<<" "<<backbone[1]<<" "<<backbone[2]<<" "<<backbone[3]<<endl;
               return 0;
            }
            for(j=0; j<4; ++j)
               backbone[j]=0;
            if( strcmp(rac[i].atom, "N")==0) backbone[0]=1;
            if( strcmp(rac[i].atom, "C")==0) backbone[1]=1;
            if( strcmp(rac[i].atom, "CA")==0) backbone[2]=1;
            if( strcmp(rac[i].atom, "O")==0) backbone[3]=1;
            
            temp=rac[i].sequence;
         }
         else
         {
            return 0;
         }
         ++i;
      }
      if(backbone[0]!=1 || backbone[1]!=1 || backbone[2]!=1 || backbone[3]!=1)
      {
         cout<<"Error Format! (Lack of Backbone Atom) End"<<endl;
         return 0;
      }
      
      return 1;
}


int ProteinStructureAnalysis::ReadProteinBackbone(array2sd *BackboneAtom)
{
      int i=0, j=0;
      vector<float> coors;
      str2vectord tmp_bb;
      
      while(strlen(rac[i].resi)>0)
      {
         if(i==0)
         {
            j=rac[i].sequence;
         }
         if(j!=rac[i].sequence)
         {
            BackboneAtom->push_back(tmp_bb);
            tmp_bb.clear();
            j=rac[i].sequence;
         }
         
         if( strcmp(rac[i].atom, "N" )==0 || strcmp(rac[i].atom, "C")==0 || 
         strcmp(rac[i].atom, "CA")==0 || strcmp(rac[i].atom, "O")==0)
         {
            string atom(rac[i].atom);
            coors.push_back(rac[i].xyz[0]);
            coors.push_back(rac[i].xyz[1]);
            coors.push_back(rac[i].xyz[2]);
            tmp_bb[atom]=coors;
            coors.clear();
         }
         else if( strcmp(rac[i].atom, "CB" )==0 && CBinclude)
         {
            string atom(rac[i].atom);
            coors.push_back(rac[i].xyz[0]);
            coors.push_back(rac[i].xyz[1]);
            coors.push_back(rac[i].xyz[2]);
            tmp_bb[atom]=coors;
            coors.clear();
         }
         ++i;
      }
     
      if(tmp_bb.size()!=0)
      {
         BackboneAtom->push_back(tmp_bb);
      }
     
      
      return 1;
}


int ProteinStructureAnalysis::ReadProteinSequence(strs *rps)
{
      int i=0, j=0;
      
      while(strlen(rac[i].resi)>0)
      {
         if(i==0)
         {
            string residue(rac[i].resi);
            rps->push_back(residue);
            j=rac[i].sequence;
         }
         if(j!=rac[i].sequence)
         {
            string residue(rac[i].resi);
            rps->push_back(residue);
            j=rac[i].sequence;
         }
         ++i;
      }
      return i;
}

int ProteinStructureAnalysis::ReadProteinSequence(strs *rps, vector<int>* index)
{
      int i=0, j=0;
      
      while(strlen(rac[i].resi)>0)
      {
         if(i==0)
         {
            string residue(rac[i].resi);
            rps->push_back(residue);
            j=rac[i].sequence;
            index->push_back(j);
         }
         if(j!=rac[i].sequence)
         {
            string residue(rac[i].resi);
            rps->push_back(residue);
            j=rac[i].sequence;
            index->push_back(j);
         }
         ++i;
      }
      //index->push_back(i);
      return i;
}

void ProteinStructureAnalysis::XYZ(int i, string &atom, vector<float> &xyz)
{
      atom.assign(rac[i].atom);
      
      xyz.push_back(rac[i].xyz[0]);
      xyz.push_back(rac[i].xyz[1]);
      xyz.push_back(rac[i].xyz[2]);
      
}
void ProteinStructureAnalysis::XYZ(int i, string &atom, string &resi, vector<float> &xyz)
{
      atom.assign(rac[i].atom);
      resi.assign(rac[i].resi);
      
      xyz.push_back(rac[i].xyz[0]);
      xyz.push_back(rac[i].xyz[1]);
      xyz.push_back(rac[i].xyz[2]);
      
}
int ProteinStructureAnalysis::Seq(int i)
{
      return rac[i].sequence;
      
}
char ProteinStructureAnalysis::Chain(int i)
{
      return rac[i].chain;
      
}
int ProteinStructureAnalysis::Atom_Num()
{
      return atom_num;
}

int ProteinStructureAnalysis::TransferData(const strs *residue_name, array2sd *cartesian)
{
      racfree();
      char ch_id='A';
      int i=0, j=0;
      str2vectord::iterator pos2;

      for(i=0; i<Max_Atom_Number_BS; i++)
      {
         rac[i].sequence=0;
         rac[i].chain=' ';
         initialization(rac[i].resi, 0, 6);
         initialization(rac[i].atom, 0, 7);
         rac[i].xyz[0]=0; rac[i].xyz[1]=0; rac[i].xyz[2]=0;
      }
      
      for(i=0, j=0; i<cartesian->size(); ++i)
      {
          //If (cartesian->at(i).size()==0) it means there is a gap.Another chain
          if(cartesian->at(i).size()==0) 
          {
             ch_id+=1;
             continue;
          }
          for(pos2=cartesian->at(i).begin(); pos2!=cartesian->at(i).end(); ++pos2)
           {
             //if(pos2->first[0]=='-' || pos2->first[0]=='+') continue;
             //if(pos2->first[0]=='H') continue;
             
             residue_name->at(i).copy(rac[j].resi, residue_name->at(i).size());
             pos2->first.copy(rac[j].atom, pos2->first.size());
             rac[j].sequence=i;
             rac[j].xyz[0]=(pos2->second).at(0);
             rac[j].xyz[1]=(pos2->second).at(1);
             rac[j].xyz[2]=(pos2->second).at(2); 
             rac[j].chain=ch_id;
             //cout<<"TD"<<rac[j].sequence<<" "<<rac[j].resi<<" "
             //<<residue_name->at(i).size()<<" "<<rac[j].atom<<" "
             //<<rac[j].xyz[0]<<" "<<rac[j].xyz[1]<<" "<<rac[j].xyz[2]<<"\n";
           
             ++j;
            }
         
       } 
       //cout<<"TransferData N="<<j<<endl;
       atom_num=j;
      /*Verify the data.
      cout<<"Verify the data:"<<endl;
      for(i=0; i<j; i++)
      {
         cout<<rac[i].sequence<<" "<<rac[j].chain<<" "<<rac[i].resi<<" "
         <<rac[i].atom<<" "<<rac[i].xyz[0]<<" "<<rac[i].xyz[1]<<" "<<rac[i].xyz[2]<<"\n";
      }
      //*/
      return atom_num;
}

int ProteinStructureAnalysis::Coordinates(array2sd *cartesian)
{
      int i=0, j=0;
      str2vectord::iterator pos2;
      
      for(i=0, j=0; i<cartesian->size(); ++i)
      {
          if(cartesian->at(i).size()==0) continue;
          for(pos2=cartesian->at(i).begin(); pos2!=cartesian->at(i).end(); ++pos2)
           {
             if(j==atom_num)
             {
                cout<<"Error in Coordinates()! Data is not matched "<<endl;
                exit(0);
             }
             
             //if(pos2->first[0]=='-' || pos2->first[0]=='+') continue;
             //if(pos2->first[0]=='H') continue;
             
             rac[j].xyz[0]=(pos2->second).at(0);
             rac[j].xyz[1]=(pos2->second).at(1);
             rac[j].xyz[2]=(pos2->second).at(2);

             //cout<<rac[j].sequence<<" "<<rac[j].resi<<" "
             //<<rac[j].atom<<" "<<rac[j].xyz[0]<<" "<<rac[j].xyz[1]<<" "<<rac[j].xyz[2]<<"\n";
           
             ++j;
            }
       } 
      
      return j;
}
/*
int ProteinStructureAnalysis::RotamerType(vector<int> &rotamer)
{
      int i=0;
      vector<int> rotemp;
      strs peptide, atoms;
      matrix resicoor;
      array2sd cartesian;
      ExtractPDB_Gapless(&peptide, &cartesian);
      for(i=0; i<peptide.size(); ++i)
      {
         if(Sockets(peptide[i], cartesian[i], resicoor, atoms)==1)
         {
            CalRotamer(peptide[i], resicoor, rotemp);
            rotamer.insert(rotamer.end(), rotemp.begin(), rotemp.end());
         }
         else
         {
            cout<<"No Rotamer Type found!\n";
            rotamer.push_back(0); rotamer.push_back(0);
            rotamer.push_back(0); rotamer.push_back(0);
         }
         //rotamer.push_back(IndexTransform(rotemp));
         //for(int j=0; j<rotemp.size(); ++j)
         //   cout<<rotemp[j]<<" ";
         //cout<<endl;
      }
      return 1;;
}
*/
int ProteinStructureAnalysis::GetDecoyData(char *input, char* file2)
{
      racfree();
      
      int i=0, j=0;
      char line[100], temp[8];
      char *x1;
      
      ifstream oi (input, ios::in);
      if(!oi)
      {
         cout<<"Error in file reading1 "<<input<<endl;
         return 0;
      }
      
      oi.seekg(0,ios::beg);
      while(oi.good())
      {
         oi.getline(line, 100, '\n');   
         x1=strstr(line, "ATOM");
         if(x1 != NULL && line[0]=='A')
         {
            copy(rac[i].resi, line, 17, 20);
            copy(rac[i].atom, line, 13, 16);
            
            initialization(temp, 0, 8);
            copy(temp, line, 12, 13);
            if(strcmp(temp, " ")!=0)
            strcat(rac[i].atom, temp);
            
            if(strcmp(rac[i].atom, "H")==0) 
            {
               //cout<<"H:"<<rac[i].atom<<" "<<rac[i].resi<<"\n";
               strcat(rac[i].atom, "N");
            }
            if(strcmp(rac[i].atom, "HG")==0 && strcmp(rac[i].resi, "SER")==0)
            {
               //cout<<"HG:"<<rac[i].atom<<" "<<rac[i].resi<<"\n";
               strcat(rac[i].atom, "1");
            }
            if(strcmp(rac[i].atom, "CD1")==0 && strcmp(rac[i].resi, "ILE")==0)
            {
               rac[i].atom[2]= '\0';
            }
            if(strcmp(rac[i].resi, "HIS")==0)
            {
               initialization(rac[i].resi, 0, 6);
               strcat(rac[i].resi, "HSD");
            }
            if(strcmp(rac[i].resi, "ILE")==0)
            {
               if(strcmp(rac[i].atom, "HD11")==0)
               {
                  rac[i].atom[3]='\0';
               }
               else if(strcmp(rac[i].atom, "HD12")==0)
               {
                  rac[i].atom[3]='\0';
                  rac[i].atom[2]='2';
               }
               else if(strcmp(rac[i].atom, "HD13")==0 )
               {
                  rac[i].atom[3]='\0';
                  rac[i].atom[2]='3';
               }
            }
            
            rac[i].chain=line[21];
            
            initialization(temp, 0, 8);
            copy(temp, line, 22, 26);
            rac[i].sequence=atoi(temp);
            /*
            initialization(temp, 0, 8);
            copy(temp, line, 30, 38);
            rac[i].xyz[0]=atof(temp);
            
            initialization(temp, 0, 8);
            copy(temp, line, 38, 46);
            rac[i].xyz[1]=atof(temp);
            
            initialization(temp, 0, 8);
            copy(temp, line, 46, 54);
            rac[i].xyz[2]=atof(temp);
            */
            i++;
            if( i>Max_Atom_Number_BS )
            {
               cout<<"Error in GetData(). Overflow from the Max_Atom_Number_BS limit!\n";
               break;
            }
         }
      }
      
      ifstream ff (file2, ios::in);
      if(!ff) 
      {
         cout<<"Error in file reading2 "<<file2<<endl;
         return 0;
      }
      int markline=0;
      float x=0, y=0, z=0;
      j=0; string tmp;
      while(ff.good())
      {
         /*
         if(markline!=0)
         {
            ff>>x;
            ff>>y;
            ff>>z;
            j++;
            cout<<x<<" "<<y<<" "<<z<<endl;
            if(j>i)
            {
               cout<<"Error in march pdb file"<<endl;
               exit(0);
            }
         }
         else
         {
            ff>>tmp;
            cout<<tmp<<" ";
         }
         
         if(tmp=="5")markline++;
         */
         if(markline==0)
         {
            ff.getline(line, 100, '\n');  
         }
         else
         {
            ff>>rac[j].xyz[0];
            ff>>rac[j].xyz[1];
            ff>>rac[j].xyz[2];
            //cout<<rac[j].xyz[0]<<" "<<rac[j].xyz[1]<<" "<<rac[j].xyz[2]<<endl;
            j++;
            if(j==i)
            {
               //cout<<"Error in march pdb file"<<endl;
               break;
            }
         }
         markline++;
      }
      if(j!=i)
      {
         cout<<"Warning "<<i<<" != "<<j<<endl;
      }
      //cout<<"Atom Number "<<i<<endl;
      /*Verify the data.
      for(j=0; j<i; j++)
      {
         cout<<"V "<<rac[j].sequence<<" "<<rac[j].resi<<" "
         <<rac[j].atom<<" "<<rac[j].xyz[0]<<" "<<rac[j].xyz[1]<<" "<<rac[j].xyz[2]<<"\n";
      }
      //*/
      x1=NULL;
      atom_num=i;
      return i;
}

int ProteinStructureAnalysis::RemoveHydrogen(array2sd &cartesian)
{
      int i=0, j=0;
      str2vectord::iterator pos2;
      for(i=0, j=0; i<cartesian.size(); ++i)
      {
         if(cartesian.at(i).size()==0) continue;
         for(pos2=cartesian.at(i).begin(); pos2!=cartesian.at(i).end(); )
         {
            if(pos2->first[0]=='-' || pos2->first[0]=='+')
            {
               cout<<"Error atom ProteinStructureAnalysis::RemoveHydrogen()"<<endl;
            }
            //cout<<i<<" "<<pos2->first<<endl;
            //Take care of using erase!!! 
            //Directly cartesian.at(i).erase(pos) will be wrong!
            if(pos2->first[0]=='H')//&&pos2->first!="HN")
            {
               //cout<<pos2->first<<endl;
               cartesian.at(i).erase(pos2++);
               ++j;
            }
            else
            {
               ++pos2;
            }

         }
       } 

      return j;
}
int ProteinStructureAnalysis::RemoveHydrogen()
{
      int i,j=0;
      for(i=0, j=0; j<atom_num; ++j)
      {
         if(rac[j].atom[0]!='H')
         {
            rac[i].sequence=rac[j].sequence;
            rac[i].chain=rac[j].chain;
            copy(rac[i].resi, rac[j].resi, 0, 6);
            copy(rac[i].atom, rac[j].atom, 0, 7);
            rac[i].xyz[0]=rac[j].xyz[0]; 
            rac[i].xyz[1]=rac[j].xyz[1]; 
            rac[i].xyz[2]=rac[j].xyz[2];
            ++i;
         }
      }
      for(j=i; j<atom_num; ++j)
      {
         rac[j].sequence=0;
         rac[j].chain=' ';
         initialization(rac[j].resi, 0, 6);
         initialization(rac[j].atom, 0, 7);
         rac[j].xyz[0]=0; rac[j].xyz[1]=0; rac[j].xyz[2]=0;
      }
      
      atom_num=i;

      return atom_num;
}

int ProteinStructureAnalysis::SubstituteResidue(const vector<int> &index, const strs &substitution)
{
      int i=0, j=0;
      for(i=0; i<atom_num; ++i)
      {
         for(j=0; j<index.size(); ++j)
         {
            if(rac[i].sequence==index[j])
            {
               if(substitution[j].size()!=3)
               {
                  cout<<"Error in ProteinStructureAnalysis::SubstituteResidue\n";
                  exit(0);
               }
               if(substitution[j]=="HIS")
               {
                  initialization(rac[i].resi, 0, 6);
                  strcat(rac[i].resi, "HSD");
                  //rac[i].resi[0]='H';
                  //rac[i].resi[1]='S';
                  //rac[i].resi[2]='D';
               }
               else
               {
                  rac[i].resi[0]=substitution[j][0];
                  rac[i].resi[1]=substitution[j][1];
                  rac[i].resi[2]=substitution[j][2];
               }
               rac[i].resi[3]='\0';
               rac[i].resi[4]='\0';
               rac[i].resi[5]='\0';
               cout<<"SubRes "<<index[j]<<" "<<rac[i].resi<<endl;
            }
         }
      }
      return 1;
}
int ProteinStructureAnalysis::RemoveSideChain()
{
      int i,j=0;
      for(i=0, j=0; j<atom_num; ++j)
      {
         if(strcmp(rac[j].atom, "C")==0 || strcmp(rac[j].atom, "CA")==0
        || strcmp(rac[j].atom, "O")==0 || strcmp(rac[j].atom, "N")==0)
         {
            rac[i].sequence=rac[j].sequence;
            rac[i].chain=rac[j].chain;
            copy(rac[i].resi, rac[j].resi, 0, 6);
            copy(rac[i].atom, rac[j].atom, 0, 7);
            rac[i].xyz[0]=rac[j].xyz[0]; 
            rac[i].xyz[1]=rac[j].xyz[1]; 
            rac[i].xyz[2]=rac[j].xyz[2];
            ++i;
         }
         
         if(strcmp(rac[j].atom, "CB")==0 && CBinclude)
         {
            rac[i].sequence=rac[j].sequence;
            rac[i].chain=rac[j].chain;
            copy(rac[i].resi, rac[j].resi, 0, 6);
            copy(rac[i].atom, rac[j].atom, 0, 7);
            rac[i].xyz[0]=rac[j].xyz[0]; 
            rac[i].xyz[1]=rac[j].xyz[1]; 
            rac[i].xyz[2]=rac[j].xyz[2];
            ++i;
         }
      }
      for(j=i; j<atom_num; ++j)
      {
         rac[j].sequence=0;
         rac[j].chain=' ';
         initialization(rac[j].resi, 0, 6);
         initialization(rac[j].atom, 0, 7);
         rac[j].xyz[0]=0; rac[j].xyz[1]=0; rac[j].xyz[2]=0;
      }
      
      atom_num=i;

      return atom_num;
}



#endif
