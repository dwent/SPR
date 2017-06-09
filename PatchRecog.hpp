//==========================================================================
// Description:
//	Interface-like patches' recognition.
//	Mainly focus on the recognition scoring function.
//
// Contact: went.dai@gmail.com
// Date   : 2015-10-20
//==========================================================================

#ifndef PATCHRECOG_H
#define PATCHRECOG_H

#include <vector>
#include <iostream>
#include <string>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <ctime>
//#include "/tjjiang/wuaiping/DomainInteraction/src/PatchRecognition/SCentropy/src/ec.hpp"


using namespace std;


//== Interface-like Patches Recognition Score Function.
//double PatchScore(vector<double> domsol,
//		     vector<vector<double > > freq,
//		     vector<vector<double > > prof,
//		     string sequence, 
//		     vector<int> patch, 
//		     vector<double> asa);
double SeqScore(vector<int> patch, 
		     string sequence, 
		     vector<double> asa);
double HydroScore(vector<int> patch, 
		     string sequence, 
		     vector<double> asa);
double ConserveScore(vector<vector<double > > freq,
		     vector<vector<double > > prof, 
		     string sequence, 
		     vector<int> patch);
double SolScore(vector<int> patch, 
		     vector<double> domsol);


//== Functions called by Score-Functions.
int ReadASAcoord(vector<double>& xca, 
		    vector<double>& yca, 
		    vector<double>& zca, 
		    string name);
int ScContactMatrix(vector<vector<double > >& SCMTX, 
		    string name);
void GenCAmatrix(vector<vector<double > >& DMTX, 
		    vector<vector<double > >& SCMTX, 
		    string asafile);
double InitPatch(vector<int>& patch, 
		    int pos, 
		    const vector<int>& si, 
		    const vector<double>& asa, 
		    const vector<vector<double > >& SCMTX);
void PatchMerge(vector<vector<int > >& allpatch, 
		    double domainasa);
int SameNum(vector<int> patch1, 
	          vector<int> patch2);


//== I/O functions and others.
int ReadID(vector<string>& allid, 
		    string IDfile);
void SubstractSS(const string sequence, 
		    const string dsspFile, 
		    string &second);
int ReadPSSM(vector<vector<double > >& freq, 
		    vector<vector<double > >& prof, 
		    string sequence, 
		    string psifile);
void MultiFace(vector<vector<string > >& multifaceID, 
		    vector<string> singleID);
bool ReadInterface(string& sequence, 
		    string& faceinfo, 
		    vector<double>& asa, 
		    string iresfile);
int ReadSol(vector<double>& domsol, 
		    string solfile);
void CompPatchFace(vector<int>& Truenum, 
		    vector<int>& Allnum, 
		    vector<int>& Prenum, 
		    vector<vector<int > > allpatch, 
		    vector<vector<int > > allface);
void OutputCovAcc(vector<vector<double > > allscore, 
		    vector<vector<int > > alltnum, 
		    vector<vector<int > > allrnum, 
		    vector<vector<int > > allpnum, 
		    vector<vector<vector<int > > > patchinfo,
		    vector<vector<double > > otherscore);
void ReadSequence(const char* fastafile,string& sequence);
int ReadRsa(string path, vector<double>& rsa);

//== Functions to be checked!
double EntropyScore(vector<int> patch, 
			string pdbfile);




#endif
