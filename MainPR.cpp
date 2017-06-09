//==========================================================================
//
// 测试使用Patch Recognition方法中的识别评分函数，并确定各项之间的权重。
//
// Contact: went.dai@gmail.com
// Date   : 2015-10-20
//
//==========================================================================


#include "PatchRecog.hpp"

using namespace std;




//======================================================================
//
// Compile: g++ -O2 -o patchrecog MainPR.cpp PatchRecog.cpp 
// Usage  : ./patchrecog ./ID ./asa ./ires ./dssp ./fasta.psi ./sol
//
//======================================================================
int main(int argc, char* argv[])
{
	
	int i=0, j=0, k=0, l=0;
	
	
	//== Read all complexes' ID
	string idfile=argv[1];
	vector<string> allcid;
	ReadID(allcid, idfile);

	
	//== Deal with the situation that one domain has multi-interfaces
	//vector<vector<string > > domfaceID;
	//MultiFace(domfaceID, allcid);
	
	
	//== read *.ires file to get sequence, interface and ASA value;
	string asapath=argv[2];
	string irespath=argv[3];
	string dssppath=argv[4];
	string psipath=argv[5];
	string solpath=argv[6];

	
	string second="";
	string sequence="";
	//string faceinfo="";
	vector<double> asa;


      clock_t start,finish;
      double totaltime;
      start=clock();

	
	double w1=1.1, w2=0.3, w3=0.1;
	
	vector<vector<double > > allscore;
	vector<vector<double > > otherscore;
	//vector<vector<int > > alltnum;
	//vector<vector<int > > allrnum;
	//vector<vector<int > > allpnum;
	//vector<vector<vector<int > > > patchinfo;
	//string stemp="";
	double tempasa=0.0;
	for(i=0; i<allcid.size(); i++)
	{
		sequence="";
		//faceinfo="";
		asa.clear();

		
		//==read fasta files of this domain
		string iresfile=irespath+allcid[i]+".fasta";
		//cout<<iresfile<<endl;
		char fastafile[100];
		strcpy(fastafile, iresfile.c_str());
		ReadSequence(fastafile, sequence);

		//sequence="GVQVETISPGDGRTFPKRGQTCVVHYTGMLEDGKKFDSSRDRNKPFKFMLGKQEVIRGWEEGVAQMSVGQRAKLTISPDYAYGATGHPGIIPPHATLVFDVELLKLE";

		//==Read rsafile
		string rsafile=asapath+allcid[i]+".rsa";
		ReadRsa(rsafile, asa);
		

/*		
		//== read all interface files (*.ires) of this domain
		string iresfile=irespath+domfaceID[i][0];
		ReadInterface(sequence, faceinfo, asa, iresfile);
		cout<<sequence<<endl;
	
		vector<string> allfaceinfo;
		allfaceinfo.push_back(faceinfo);
		if(domfaceID[i].size()>1)
		{
			for(j=1; j<domfaceID[i].size(); j++)
			{
				string tempseq="";
				string tempface="";
				vector<double> tempasa;
				string tempfile=irespath+domfaceID[i][j];
				bool faceflag=ReadInterface(tempseq, tempface, tempasa, tempfile);
				if(!faceflag) continue;
				else allfaceinfo.push_back(tempface);
			}
		}
		if(allfaceinfo.empty()) continue;
*/		


		
		//== read secondary structure information
		second="";
		string dsspfile=dssppath+allcid[i]+".dssp";
		SubstractSS(sequence, dsspfile, second);
		if(second.size()!=sequence.size()) continue;
		
		
		//== read sequence profile information
		string psifile=psipath+allcid[i]+".fasta.psi";
		vector<vector<double > > freq;
		vector<vector<double > > prof;
		ReadPSSM(freq, prof, sequence, psifile);
		if(freq.size()!=sequence.size()) continue;


		//== solvation energy.
		string solfile=solpath+allcid[i]+".sol";
		vector<double> domsol;
		ReadSol(domsol, solfile);
		if(domsol.size()!=sequence.size()) continue;
					

/*
		//== Recognize all the interface(s) positions in the same domain
		vector<vector<int > > allface;
		vector<int> itemp;
		for(j=0; j<allfaceinfo.size(); j++)
		{
			itemp.clear();
			for(k=0; k<allfaceinfo[j].size(); k++)
			{
				if(allfaceinfo[j].substr(k,1)=="1")
				{
					itemp.push_back(k);
				}
			}
			
			if(itemp.size()>9) allface.push_back(itemp);
		}
		if(allface.empty()) continue;
*/
		

		//== read *.asa file to get XYZ-coordinates of all CA atoms
		string asafile=asapath+allcid[i]+".asa";
		int seqnum=asa.size();
		vector<vector<double > > DMTX;
		vector<vector<double > > SCMTX;
		DMTX.resize(seqnum);
		SCMTX.resize(seqnum);
		for(j=0; j<seqnum; j++)
		{
			DMTX[j].resize(seqnum);
			SCMTX[j].resize(seqnum);
		}
		GenCAmatrix(DMTX, SCMTX, asafile);

		
		//== identify all surface residues in the domain structure
		vector<int> si;
		for(j=0; j<asa.size(); j++)
		{
			if(asa[j]>1.0) si.push_back(j);
		}
		
		
		//== generate patches to simulate real interfaces
		vector<vector<int > > allpatch;
		for(j=0; j<si.size(); j++) //generate si.size() number of patches.
		{
			vector<int> tpatch;
			tempasa=InitPatch(tpatch, j, si, asa, SCMTX);
				
			sort(tpatch.begin(), tpatch.end());
			allpatch.push_back(tpatch);
		}
		if(allpatch.empty()) continue;
		
		double domainasa=0.0;
		for(j=0; j<asa.size(); j++) domainasa += asa[j];
		PatchMerge(allpatch, domainasa);
		
/*			
		//== calculate the maximum residues coverage and precision of each patch to all interface(s).
		vector<int> truenum, resinum, prenum;
		CompPatchFace(truenum, resinum, prenum, allpatch, allface);
		alltnum.push_back(truenum);
		allrnum.push_back(resinum);
		allpnum.push_back(prenum);
		patchinfo.push_back(allpatch);
*/
		
		
		//== Calculate the patch score
		vector<double> score;
		for(j=0; j<allpatch.size(); j++)
		{
             	double Eseq=SeqScore(allpatch[j], sequence, asa);
            	double Ehydro=HydroScore(allpatch[j], sequence, asa);
            	double Econs=ConserveScore(freq, prof, sequence, allpatch[j]);
            	double Esol=SolScore(allpatch[j], domsol);
			double tempscore=Eseq+w1*Ehydro+w2*Econs+w3*Esol;
			score.push_back(tempscore);
		}
		allscore.push_back(score);
		
		cout<<allcid[i]<<endl;
		for(j=0; j<allpatch.size(); j++ )
		{
			cout<<"Patch"<<j<<": "<<score[j]<<endl;
			for(k=0; k<allpatch[j].size(); k++)
			{
				cout<<allpatch[j][k]<<'\t';
			}cout<<endl;
		}
	}

      finish=clock();
      totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
      cout<<"\nTime is"<<totaltime<<"second！"<<endl;


	
	return 0;
}


