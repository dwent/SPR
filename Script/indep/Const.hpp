#ifndef Consthpp
#define Consthpp
# include "Directory.hpp"
# include <map>
# include <string>
# include <vector>
#include <cstdlib>
#include <cstring> 
using namespace std;

const float PI=3.1415926535897932384626433832795028841971693993751;
const float PIx2=2*PI;
const float Extra=1.0e-4;
const float UpMax=1.0e+10;

typedef multimap<string, vector <float> > Rotlib;
typedef multimap<string, int> strs2int;
typedef map<string, vector <string> > str2strs;
typedef map<string, vector <float> > str2vectord;
typedef map<string, int> str2int;
typedef map<string, float> str2float;
typedef vector <str2vectord> array2sd;
typedef vector <string> strs;
typedef vector<vector<float> > matrix;
typedef vector<vector<int> > matxint;
typedef vector<float> myvector;
typedef vector<int> myvecint;

struct ResiAtomCoor_BS
{
      int sequence;
      char chain;
      char resi[6];
      char atom[7];
      float xyz[3];
};

enum Residues 
      {
      GLY, ALA, VAL, LEU, 
      ILE, MET, PRO, CYS, 
      PHE, TYR, TRP, ARG, 
      LYS, HIS, ASP, GLU, 
      SER, THR, ASN, GLN
      };
//#define CB     // CB include in main chain or not.
//#define AMBER  // Using AMBER data. If not, Using CHARMM data.


//#define  LIBP /tjjiang/caoyang/HA/2008-11-5/basiclib
//const string pan_para="/tjjiang/caoyang/Cheer/HA/2008-11-5/scoring/implict-solv/pan/sphere2.raw";
//const string top_para="/tjjiang/caoyang/Cheer/HA/2008-11-5/basiclib/top_rotamer_amber_h.inp";
//const string rot_para="/tjjiang/caoyang/Cheer/HA/2008-11-5/RotamerLibrary/rot.txt";
//const string bbdeplib="/tjjiang/caoyang/Cheer/HA/2008-11-5/RotamerLibrary/bbdep02.May.lib";
//const string ptcharmm="/tjjiang/caoyang/Cheer/HA/2008-11-5/minicharmm/";
#endif





/*
/////////////////NoUse///////////////////
int readconfig(const char *configfile)
{
      string mark;
      ifstream FILE(configfile, ios::in);
      if(!FILE) {cout<<"Error config File "<<configfile<<endl; exit(0);}
      cout<<"Reading config file "<<configfile<<endl;
      while(FILE>>mark)
      {
         if(mark=="pan") FILE>>pan_para;
         else if(mark=="top") FILE>>top_para;
         else if(mark=="rot") FILE>>rot_para;
         else if(mark=="lib") FILE>>bbdeplib;
         else cout<<"Text "<<mark<<endl;
      }
      FILE.close();
}
/////////////////////////////////////////
*/









