//==========================================================================
// Description:
//	Interface-like patches' recognition.
//	Mainly focus on the recognition scoring function.
//
//Modify:ReadSol() line983, str2!="*",eg. in different chains will have "*" in file
//Modify2:SeqScore() is modified by the new SeqDistr from SCOPe2.05
//
// Contact: went.dai@gmail.com
// Date   : 2015-05-28
//==========================================================================


#include "PatchRecog.hpp"

using namespace std;



//=============================================================================
// Calculate the score of the specific Domain_Surface_Patch;
// Formula:
//		E = E(seq) + a*E(hydro) + b*E(conserved).
//=============================================================================
//double PatchScore(vector<double> domsol,
//		     vector<vector<double > > freq,
//		     vector<vector<double > > prof,
//		     string sequence, 
//		     vector<int> patch, 
//		     vector<double> asa)
//{
	
//	double Eseq=SeqScore(patch, sequence, asa);
//	double Ehydro=HydroScore(patch, sequence, asa);
//	double Econs=ConserveScore(freq, prof, sequence, patch);
//	double Esol=SolScore(patch, domsol);
	
	//double score=Eseq+1.1*Ehydro+0.4*Econs;	//local optimal.
	//double score=Eseq+1.4*Ehydro+0.35*Econs;	//refined weights.
//	double score=Eseq+1.4*Ehydro+0.35*Econs+0.25*Esol;
	
//	return score;
//}


//===================================================================
// Calculate the score of patch by its sequence bias
//===================================================================
double SeqScore(vector<int> patch, 
		     string sequence, 
		     vector<double> asa)
{
	int i=0, j=0;

	char RNAME[21]="ARNDCQEHILKMFPSTWYVG";
	double SeqDistr[20]={0.841,1.346,0.958,0.830,1.172,0.909, 0.805,1.147,1.084,1.144,0.784,1.451,1.334,1.109,0.873,0.966,1.284,1.368,0.994,0.823};
	double RefASA[20]={27.8,94.7,60.1,60.6,15.5,68.7,68.2,50.7,22.8,27.6,103.0,33.5,25.5,51.5,42.0,45.0,34.7,55.2,23.7,24.5}; //H JANJ780101 //D Average accessible surface area (Janin et al., 1978)		

	const char* chseq=sequence.c_str();		
	double score=0.0;
	double allasa=0.0;
	for(i=0; i<patch.size(); i++)
	{
		for(j=0; j<20; j++)
		{
			if(chseq[patch[i]]==RNAME[j])
			{
				score += SeqDistr[j]*asa[patch[i]]/RefASA[j];
				break;
			}
		}
		
	}

	return score;
}


//===========================================================================
// Calculate the score of patch by its residue hydrophobic bias
//===========================================================================
double HydroScore(vector<int> patch, 
		     string sequence, 
		     vector<double> asa)
{
	int i=0, j=0;

	char RNAME[21]="ARNDCQEHILKMFPSTWYVG";
	double Hydro[20]={0.2,-0.7,-0.5,-1.4,1.9,-1.1,-1.3,0.4,1.4,0.5,-1.6,0.5,1.0,-1.0,-0.7,-0.4,1.6,0.5,0.7,-0.1}; //CASG920101 Hydrophobicity scale from native protein structures (Casari-Sippl, 1992)
	double RefASA[20]={27.8,94.7,60.1,60.6,15.5,68.7,68.2,50.7,22.8,27.6,103.0,33.5,25.5,51.5,42.0,45.0,34.7,55.2,23.7,24.5}; //H JANJ780101 //D Average accessible surface area (Janin et al., 1978)		

	const char* chseq=sequence.c_str();		
	double score=0.0;
	double allasa=0.0;
	for(i=0; i<patch.size(); i++)
	{
		for(j=0; j<20; j++)
		{
			if(chseq[patch[i]]==RNAME[j])
			{
				score += Hydro[j];
				break;
			}
		}
		
	}

	return score;
}


//===========================================================================
// Calculate the score of a patch by its sequence conservation
//===========================================================================
double ConserveScore(vector<vector<double > > freq,
		     vector<vector<double > > prof, 
		     string sequence, 
		     vector<int> patch)
{
	double score=0.0;
	int i=0, j=0, k=0;
	
	char RNAME[21]="ARNDCQEGHILKMFPSTWYV";
	double BLOSUM62[20]={4,5,6,6,9,5,5,6,8,4,4,5,5,6,7,4,5,11,7,4};
	
	const char* chseq=sequence.c_str();
	double tempscore=0.0;
	for(i=0; i<patch.size(); i++)
	{
		int ti=patch[i];
		
		tempscore=0.0;
		//for(j=0; j<20; j++)	//第一种保守性计算方法（这里无效）；
		//{
		//	if(freq[ti][j]>1.0e-10) tempscore += (freq[ti][j]/100.0)*log(freq[ti][j]/100.0);
		//}
		
		for(j=0; j<20; j++)
		{
			if(chseq[ti]==RNAME[j])
			{
				//if(prof[ti][j]-BLOSUM62[j]>0) tempscore=prof[ti][j]-BLOSUM62[j];	//第二种保守性计算方法（这里无效）；
				//else tempscore=0.0;
				tempscore=prof[ti][j]-BLOSUM62[j];	//第三种保守性计算方法（有效的方法）；
			
				break;
			}
		}
		
		score += tempscore;
	}
	
	return score;
}


//===========================================================================
// Calculate the score of a patch by its solvation energy.
//===========================================================================
double SolScore(vector<int> patch, 
		     vector<double> domsol)
{
	int i=0;

	double score=0.0;
	for(i=0; i<patch.size(); i++)
	{
		score += domsol[patch[i]];
	}

	return score;	
}



//######################################################################


//===================================================================
// Read all XYZ-coordinates of CA atoms in *.asa file
//===================================================================
int ReadASAcoord(vector<double>& xca, 
		    vector<double>& yca, 
		    vector<double>& zca, 
		    string name)
{
	int i=0, j=0, k=0;
	
	ifstream fin(name.c_str());
	if( !fin.is_open() )
	{
		cout << "Error: can not open " << name << endl;
		return -2;
	}
	
	
	char buf[200];
	while( !fin.eof() )
	{
		fin.getline(buf, 200);
		string head="";
		head.push_back(buf[0]);
		head.push_back(buf[1]);
		head.push_back(buf[2]);
		head.push_back(buf[3]);
		
		if(head=="ATOM")
		{
			string coord="";
			string smark="";
			smark.push_back(buf[13]);
			smark.push_back(buf[14]);
			
			if(smark=="CA")
			{
				coord="";
				for(i=30; i<39; i++)
				{
					if(buf[i]!=' ') coord.push_back(buf[i]);
				}
				xca.push_back( atof(coord.c_str()) );
				
				coord="";
				for(i=39; i<47; i++)
				{
					if(buf[i]!=' ') coord.push_back(buf[i]);
				}
				yca.push_back( atof(coord.c_str()) );				

				coord="";
				for(i=47; i<55; i++)
				{
					if(buf[i]!=' ') coord.push_back(buf[i]);
				}
				zca.push_back( atof(coord.c_str()) );				
			}
		}
	}
	fin.close();

	
	return 0;
}


//=======================================================================
// Get the side-chain contact matrix of all residue-paires
//=======================================================================
int ScContactMatrix(vector<vector<double > >& SCMTX, 
		    string name)
{
	int i=0, j=0, k=0, l=0;
	
	ifstream fin(name.c_str());
	if( !fin.is_open() )
	{
		cout << "Error: can not open " << name << endl;
		return -1;
	}
	
	
	vector<vector<double > > sccoord;
	vector<double> vtemp;
	
	char buf[200];
	string resimark="", oldmark="";
	int index=0;
	while( !fin.eof() )
	{
		fin.getline(buf, 200);
		string head="";
		head.push_back(buf[0]);
		head.push_back(buf[1]);
		head.push_back(buf[2]);
		head.push_back(buf[3]);
		
		resimark="";
		for(i=17; i<27; i++) resimark.push_back(buf[i]);
		
		
		if(head=="ATOM")
		{
			if(index>0 && resimark!=oldmark)
			{
				sccoord.push_back(vtemp);
				vtemp.clear();
			}
			
			string coord="";
			string smark="";
			smark.push_back(buf[13]);
			smark.push_back(buf[14]);
			
			if(smark!="N " && smark!="CA" && smark!="C " && smark!="O ")
			{
				coord="";
				for(i=30; i<39; i++)
				{
					if(buf[i]!=' ') coord.push_back(buf[i]);
				}
				vtemp.push_back( atof(coord.c_str()) );
				
				coord="";
				for(i=39; i<47; i++)
				{
					if(buf[i]!=' ') coord.push_back(buf[i]);
				}
				vtemp.push_back( atof(coord.c_str()) );				

				coord="";
				for(i=47; i<55; i++)
				{
					if(buf[i]!=' ') coord.push_back(buf[i]);
				}
				vtemp.push_back( atof(coord.c_str()) );				
			}
			
			if(resimark.substr(0,3)=="GLY" && smark=="CA")
			{
				coord="";
				for(i=30; i<39; i++)
				{
					if(buf[i]!=' ') coord.push_back(buf[i]);
				}
				vtemp.push_back( atof(coord.c_str()) );
				
				coord="";
				for(i=39; i<47; i++)
				{
					if(buf[i]!=' ') coord.push_back(buf[i]);
				}
				vtemp.push_back( atof(coord.c_str()) );				

				coord="";
				for(i=47; i<55; i++)
				{
					if(buf[i]!=' ') coord.push_back(buf[i]);
				}
				vtemp.push_back( atof(coord.c_str()) );				
				
			}
		}
		
		oldmark=resimark;
		index++;
	}
	fin.close();
	sccoord.push_back(vtemp);
	
	
	
	//=================================
	double cutoff=25.0*25.0;
	double dist=0.0;
	double x1=0.0, y1=0.0, z1=0.0;
	double x2=0.0, y2=0.0, z2=0.0;
	for(i=0; i<sccoord.size(); i++)
	{
		for(j=i+1; j<sccoord.size(); j++)
		{
			bool flag=false;
			double minD=100.0*100.0;
			for(k=0; k<sccoord[i].size()/3; k++)
			{
				x1=sccoord[i][3*k+0];
				y1=sccoord[i][3*k+1];
				z1=sccoord[i][3*k+2];
				
				for(l=0; l<sccoord[j].size()/3; l++)
				{
					x2=sccoord[j][3*l+0];
					y2=sccoord[j][3*l+1];
					z2=sccoord[j][3*l+2];
					
					dist = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
					if(dist<minD) minD=dist;
					
					if(dist>cutoff)
					{
						SCMTX[i][j]=sqrt(cutoff);
						SCMTX[j][i]=sqrt(cutoff);
					
						flag=true;
						break;
					}
				}
				
				if(flag) break;
				else
				{
					SCMTX[i][j]=sqrt(minD);
					SCMTX[j][i]=sqrt(minD);
				}
			}
		}
	}
	
	
	
	return 0;
}


//=======================================================================
// Generate the N*N distance matrix of all CA-CA atom-paires
//=======================================================================
void GenCAmatrix(vector<vector<double > >& DMTX, 
		    vector<vector<double > >& SCMTX, 
		    string asafile)
{
	int j=0, k=0;
	
	vector<double> xca;
	vector<double> yca;
	vector<double> zca;	
	ReadASAcoord(xca, yca, zca, asafile);

	int seqnum=DMTX.size();
	for(j=0; j<seqnum; j++)
	{
		DMTX[j][j]=0.0;
		for(k=j+1; k<seqnum; k++)
		{
			double dist=sqrt( (xca[j]-xca[k])*(xca[j]-xca[k])+(yca[j]-yca[k])*(yca[j]-yca[k])+(zca[j]-zca[k])*(zca[j]-zca[k]) );
			DMTX[j][k]=dist;
			DMTX[k][j]=dist;
		}
	}
	
	
	ScContactMatrix(SCMTX, asafile);	
}


//===================================================================
// Generate surface patch with circle-expanding method
//===================================================================
double InitPatch(vector<int>& patch, 
		    int pos, 
		    const vector<int>& si, 
		    const vector<double>& asa, 
		    const vector<vector<double > >& SCMTX)
{
	
	int fi=si[pos];		//the original point (residue).
	patch.push_back(fi);
	
	double patchasa=0.0;
	patchasa += asa[fi];
	
	//record some shortest distance residues to point fi in surface
	for(int i=0; i<si.size(); i++)
	{
		int ti=si[i];
		if(SCMTX[fi][ti]>15.0) continue;
		else if( SCMTX[fi][ti]>2.0 && SCMTX[fi][ti]<5.0)
		{
			patch.push_back(ti);
			patchasa += asa[ti];
		}
		else if(SCMTX[fi][ti]>5.0 && SCMTX[fi][ti]<7.0 && asa[ti]>20.0)
		{
			patch.push_back(ti);
			patchasa += asa[ti];
		}
		else if(SCMTX[fi][ti]>7.0 && SCMTX[fi][ti]<9.0 && asa[ti]>40.0)
		{
			patch.push_back(ti);
			patchasa += asa[ti];
		}
		else if(SCMTX[fi][ti]>9.0 && SCMTX[fi][ti]<11.0 && asa[ti]>60.0)
		{
			patch.push_back(ti);
			patchasa += asa[ti];
		}
		else if(SCMTX[fi][ti]>11.0 && SCMTX[fi][ti]<13.0 && asa[ti]>80.0)
		{
			patch.push_back(ti);
			patchasa += asa[ti];
		}
		else if(SCMTX[fi][ti]>13.0 && SCMTX[fi][ti]<15.0 && asa[ti]>100.0)
		{
			patch.push_back(ti);
			patchasa += asa[ti];
		}
	}

	return patchasa;
}


//================================================================
// Merge patches into a integrated patch when their
//	 residues overlap large than a cutoff.
//================================================================
void PatchMerge(vector<vector<int > >& allpatch, 
		    double domainasa)
{
	
	int i=0, j=0, k=0, l=0;
	double cutoff=0.9;
	if(domainasa<5000) cutoff=0.8;
	else if(domainasa<7500) cutoff=0.7;
	else if(domainasa<10000) cutoff=0.6;
	else cutoff=0.5;
	
	
	//== cluster all patches according to their overlaped-rate
	vector<vector<int > > clusterid;
	vector<int> vtemp;
	vtemp.push_back(0);
	clusterid.push_back(vtemp);
	vtemp.clear();
	for(i=1; i<allpatch.size(); i++)
	{
		bool addflag=true;
		for(j=0; j<clusterid.size(); j++)
		{
			bool breakflag=false;
			int samenum=0;
			for(k=0; k<clusterid[j].size(); k++)
			{
				samenum=SameNum(allpatch[i], allpatch[clusterid[j][k]]);
				if(samenum*2.0/(allpatch[i].size()+allpatch[clusterid[j][k]].size())<cutoff)
				{
					breakflag=true;
					break;
				}
			}
			
			if(!breakflag)
			{
				clusterid[j].push_back(i);
				addflag=false;
				break;
			}
		}
		
		if(addflag)
		{
			vtemp.clear();
			vtemp.push_back(i);
			clusterid.push_back(vtemp);
		}
	}
	
	
	//== Merge every cluster into a integrated large patch
	vector<vector<int > > newpatch;
	vtemp.clear();
	for(i=0; i<clusterid.size(); i++)
	{
		vtemp=allpatch[clusterid[i][0]];
		for(j=1; j<clusterid[i].size(); j++)
		{
			for(k=0; k<allpatch[clusterid[i][j]].size(); k++)
			{
				bool flag=false;
				for(l=0; l<vtemp.size(); l++)
				{
					if(allpatch[clusterid[i][j]][k]==vtemp[l])
					{
						flag=true;
						break;
					}
				}
				
				if(!flag) vtemp.push_back(allpatch[clusterid[i][j]][k]);
			}
		}
		
		newpatch.push_back(vtemp);
		vtemp.clear();
	}
	

	allpatch.clear();
	allpatch=newpatch;
}


//=========================================================================
// Calculate the same number of residues in patch1 and patch2
//=========================================================================
int SameNum(vector<int> patch1, vector<int> patch2)
{
	int i=0, j=0;
	int samenum=0;
	
	for(i=0; i<patch1.size(); i++)
	{
		for(j=0; j<patch2.size(); j++)
		{
			if(patch1[i]==patch2[j])
			{
				samenum++;
				break;
			}
		}
	}
	
	return samenum;
}


//#########################################################################


//===================================================================
// Read all IDs from the ID-list file
//===================================================================
int ReadID(vector<string>& allid, string IDfile)
{
	ifstream fcid;
	fcid.open(IDfile.c_str());
	if( !fcid.is_open() )
	{
		cout << "Can not open " << IDfile << endl;
		exit(-1); 
	}


	string stemp="";
	while( !fcid.eof() )
	{
		stemp="";
		fcid >> stemp;
		if(stemp!="") allid.push_back(stemp);
	}
	fcid.close();
	
	return 0;	
}


//=============================================================================
// Extract secondary structure information(H,E,C) from DSSP file
//=============================================================================
void SubstractSS(const string sequence, 
		    const string dsspFile, 
		    string &second)
{
	string tempseq="";
	
	string stemp="";
	ifstream infileDssp(dsspFile.c_str());
	if(!infileDssp) cerr<<"error when open dssp file"<<endl;
		
	do
	{
		infileDssp>>stemp;
		if(stemp!="#") 
		{
			getline(infileDssp,stemp);
			continue;
		}
		else
		{
			getline(infileDssp,stemp);
			break;
		}
		
	}while(!infileDssp.eof());
				
	do
	{
		stemp.clear();
		getline(infileDssp,stemp);
		
		if(stemp.size()==0) break;

		if(stemp[13]=='!')
		{
			//if(second!="")
			//{
			//	const char* str=second.substr(second.size()-1,1).c_str();
			//	second.push_back(str[0]);
			//}
			//else second.push_back('C');
			continue;
		}
		else tempseq.push_back(stemp[13]);
		
		if(stemp[10]!=' ') continue;
			
		if((stemp[16]=='H')||(stemp[16]=='G')||(stemp[16]=='I'))
		{
			second.push_back('H');
		}
		else if(stemp[16]=='E')
		{
			second.push_back('E');
		}
		else
		{
			second.push_back('C');
		}
		
	}while(!infileDssp.eof());
	infileDssp.close();

	if(sequence.size()<tempseq.size())
	{
		int i=1;
		while(i<sequence.size() && sequence.size()<tempseq.size())
		{
			if(sequence.substr(i,1)!=tempseq.substr(i,1))
			{
				second.erase(second.begin()+i);
				tempseq.erase(tempseq.begin()+i);
			}
			else i++;
		}
	}

}


//======================================================================
// Read PSSM profile that generated from PSI-BLAST.
//======================================================================
int ReadPSSM(vector<vector<double > >& freq, 
		    vector<vector<double > >& prof, 
		    string sequence, 
		    string psifile)
{
	const char* pssmfile=psifile.c_str();
	int i=0, j=0, k=0;
	char buf[200];
	
	ifstream fin(pssmfile);
	if( !fin.is_open() )
	{
		cout << "\nReadPSSM()--Error: can not open " << pssmfile << endl;
		exit(-1);
	}
	
	string tempseq="";
	int index=0;
	vector<double> tempscore;
	vector<double> tempfreq;
	while( !fin.eof() )
	{
		strcpy(buf, "");
		fin.getline(buf, 200);
		
		if( index>=3 && strlen(buf)>=150 )
		{
			tempscore.clear();
			tempfreq.clear();
			
			//get sequence in *.psi file
			tempseq.push_back(buf[6]);
			
			//get score
			i=9;
			for(j=0; j<20; j++)
			{
				i=9+3*j;
				while( buf[i]==' ' )
					i++;
				tempscore.push_back( atof(&buf[i]) );
				
			}
			prof.push_back(tempscore);
			
			//get freq
			i=70;
			for(j=0; j<20; j++)
			{
				i=70+4*j;
				while( buf[i]==' ' )
					i++;
				tempfreq.push_back( atof(&buf[i])/100.0 );
			}
			freq.push_back(tempfreq);
		}
		
		index++;
	}
	fin.close();
	
	
	if(sequence.size()!=tempseq.size())	//ATTENTATION: normalize the length of structure and profile!
	{
		int pos=0;
		for(i=0; i<tempseq.size()-4; i++)
		{
			if(sequence.substr(0,4)==tempseq.substr(i,4))
			{
				pos=i;
				break;
			}
		}
		
		if(pos>0)
		{
			for(j=pos-1; j>=0; j--)
			{
				tempseq.erase(tempseq.begin()+j);
				freq.erase(freq.begin()+j);
				prof.erase(prof.begin()+j);
			}
		}
		
		//cout << sequence << endl;
		//cout << tempseq << endl << endl;
	}
	
	tempscore.clear();
	tempfreq.clear();
		
	return 0;
}


//==========================================================================
// Read sequence from fastafile;
//==========================================================================
void ReadSequence(const char* fastafile,string& sequence)
{
	ifstream fin(fastafile);
	if( !fin.is_open() )
	{
		cout << "\nCan not open fasta file in ReadFasta()!!\n";
		exit(-1);	
	}
	
	char buf[5000]; 
	char seq[5000]; 
	strcpy(seq, ""); 
	int seqnum=0;
	while( !fin.eof() )
	{
		strcpy(buf, "");
		fin.getline(buf, 5000);
		
		if(buf[0]!='>' && strlen(buf)>0)
		{
			strcat(seq, buf);
		}
	}
	
	seqnum=strlen(seq); 
	for(int i=0; i<seqnum; i++)
	{
		sequence.push_back(seq[i]);
	}
	
	fin.close(); 
}


//=================================================================
//Read rsa from rsafile
//=================================================================
int ReadRsa(string path, vector<double>& rsa)
{
	int i=0, j=0;
	
	ifstream fin;
	fin.open(path.c_str());
	if(!fin.is_open())
	{
		cout << "Error: can not open " << path << endl;
		return -1;
	}
	
	char buf[200]="";
	//string tmark="";
	while(!fin.eof())
	{
		strcpy(buf, "");
		fin.getline(buf, 200);
		if(strncmp(buf,"RES",3)==0)
		{
			//tmark="";
			//for(i=0; i<14; i++) tmark.push_back(buf[i]);
			
			//mark.push_back(tmark);
			
			j=13;
			while(buf[j]==' ') j++;
			double val=atof(&buf[j]);
			rsa.push_back(val);
		}
	}
	fin.close();
	
	
	return 0;
}


//====================================================================================
// If a domain has multi-interfaces, organize them into a same vector.
//====================================================================================
void MultiFace(vector<vector<string > >& multifaceID, 
		    vector<string> singleID)
{
	int i=0, j=0;
	
	vector<string> vtemp;
	i=0;
	while(i<singleID.size())
	{
		vtemp.clear();
		vtemp.push_back(singleID[i]);
		
		for(j=i+1; j<singleID.size(); j++)
		{
			if(singleID[j].substr(0,7)==vtemp[0].substr(0,7))
			{
				vtemp.push_back(singleID[j]);
			}
			else break;
		}
		
		multifaceID.push_back(vtemp);
		i += vtemp.size();
	}	
}


//====================================================================================
// If a domain has multi-interfaces, organize them into a same vector.
//====================================================================================
bool ReadInterface(string& sequence, 
		    string& faceinfo, 
		    vector<double>& asa, 
		    string iresfile)
{
	ifstream fid;
	fid.open(iresfile.c_str());
	if( !fid.is_open() )
	{
		cout << "Warning: can not open " << iresfile << endl;
		return false;
	}
	
	string stemp="";
	bool flag=false;
	while( !fid.eof() )
	{
		stemp="";
		
		fid >> stemp;
		if(stemp=="SEQ:") fid >> sequence;
		else if(stemp=="INT:") fid >> faceinfo;
		else if(stemp=="RSA:")
		{
			flag=true;
			continue;
		}
		
		if(flag)
		{
			if(stemp!="") asa.push_back(atof(stemp.c_str()));
		}
	}		
	fid.close();
	
	return true;		
}


//===========================================================
// Read residue solvation energy from input file.
//===========================================================
int ReadSol(vector<double>& domsol, 
		    string solfile)
{
	ifstream fin;
	fin.open(solfile.c_str());
	if( !fin.is_open() )
	{
		cout << "Can not open " << solfile << endl;
		return -1; 
	}

	string str1="", str2="", str3="", str4="";
	while( !fin.eof() )
	{
		str3="";
		fin >> str1 >> str2 >> str3 >> str4;
		if(str3!="" && str2!="*") domsol.push_back(atof(str3.c_str()));
	}
	fin.close();
	
	return 0;	

}



//===================================================================================
// Compare the maximum overlap between every patch and all interfaces
//===================================================================================
void CompPatchFace(vector<int>& Truenum, 
		    vector<int>& Allnum, 
		    vector<int>& Prenum, 
		    vector<vector<int > > allpatch, 
		    vector<vector<int > > allface)
{
	int j=0, m=0, l=0, k=0;
	int num=0;
	for(j=0; j<allpatch.size(); j++)
	{
		num=0;
		double maxsnum=-1.0e10;
		int tnum=0;
		int mi=-1;
		for(m=0; m<allface.size(); m++)
		{
			int samenum=0;
			for(l=0; l<allface[m].size(); l++)
			{
				for(k=0; k<allpatch[j].size(); k++)
				{
					if(allpatch[j][k]==allface[m][l])
					{
						samenum++;
						break;
					}
				}
			}
			
			if(samenum*samenum*1.0/(allface[m].size()*allpatch[j].size())>maxsnum)
			{
				maxsnum=samenum*samenum*1.0/(allface[m].size()*allpatch[j].size());
				tnum=samenum;
				mi=m;
			}
			
			num+=allface[m].size();
		}
		
		Truenum.push_back(tnum);
		Allnum.push_back(num);
		Prenum.push_back(allpatch[j].size());
	}
	
}


//===========================================================================
// Output the Coverage and Accuracy of the Patch Recognition
//	 method for the input testing database.
//===========================================================================
void OutputCovAcc(vector<vector<double > > allscore, 
		     vector<vector<int > > alltnum, 
		     vector<vector<int > > allrnum, 
		     vector<vector<int > > allpnum, 
		     vector<vector<vector<int > > > patchinfo,
		     vector<vector<double > > otherscore)
{
	int i=0,j=0,k=0,l=0;
	
	//for(double tn=-10.0; tn<10.0; tn=tn+0.05)
	//{
	
	vector<double> vd;
	int topN=2;
	vector<int> allopti;
	int vtnum=0, vrnum=0, vpnum=0;
	for(i=0; i<allscore.size(); i++)
	{
		vd.clear();
		vd=allscore[i];
		//for(j=0; j<allscore[i].size(); j++) vd.push_back(allscore[i][j]+tn*otherscore[i][j]);

		int mi=(int) (max_element(vd.begin(), vd.end())-vd.begin());
		
		vtnum += alltnum[i][mi];
		vrnum += allrnum[i][mi];
		vpnum += allpnum[i][mi];
		allopti.clear();
		allopti.push_back(mi);
		
		if(topN>1)
		{
			int testnum=0;
			vd[mi]=-1.0e10;
			j=1;
			while(j<topN)
			{
				bool flag=true;
				mi=(int) (max_element(vd.begin(), vd.end())-vd.begin());
				for(int n=0; n<allopti.size(); n++)
				{
					int opti=allopti[n];
					int samenum=0;
					for(k=0; k<patchinfo[i][opti].size(); k++)
					{
						for(l=0; l<patchinfo[i][mi].size(); l++)
						{
							if(patchinfo[i][opti][k]==patchinfo[i][mi][l])
							{
								samenum++;
								break;
							}
						}
						
						if(samenum>0)
						{
							flag=false;
							break;
						}
					}
					
					if(!flag) break;
				}
				
				if(!flag) vd[mi]=-1.0e10;
				else
				{
					vtnum += alltnum[i][mi];
					vpnum += allpnum[i][mi];
					allopti.push_back(mi);
					vd[mi]=-1.0e10;
					
					j++;
				}
				
				testnum++;
				if(testnum>=patchinfo[i].size()) break;
			}
		}
		
	}
	
	//cout << "\ntn= " << tn << endl;	
	cout << "Coverage: " << '\t'<<vtnum << "/" << vrnum << " = " <<vtnum*1.0/vrnum;
	cout <<'\t'<< "Precision : " << '\t'<<vtnum << "/" << vpnum << " = " << vtnum*1.0/vpnum;
      cout <<'\t'<<"weight:" <<'\t'<<(vtnum*1.0/vrnum)*(vtnum*1.0/vpnum)<< endl;
	
	//}
}


//######################################################################3


//===========================================================================
// Calculate the score of the side-chain conformational entropy
//	 of a patch.
//===========================================================================
double EntropyScore(vector<int> patch, 
			string pdbfile)
{
	double score=0.0;
	int i=0, j=0;
	const char* chfile=pdbfile.c_str();
	char infile[1000]="";
	strcpy(infile, chfile);
	
/*	
	vector<int> Inum, Rnum;
	EC bsc;			
	bsc.GetPDB(infile);
	bsc.GetSub();
	bsc.GetRot("./SCentropy/data/Rot70.dat");
	bool flag=bsc.PhiPsi();
	
	if(!flag) return -1.0e10;
		
	bsc.GetLib("./SCentropy/data/BBDepSort.bin");
	bsc.GetTop("./SCentropy/data/Miao.top");
	bsc.Build();	
	bsc.SelfEnergy(Inum, Rnum);
	
	for(i=0; i<patch.size(); i++)
	{
		j=patch[i];
		
		score += (Rnum[j]*1.0/Inum[j])*log(Rnum[j]*1.0/Inum[j]);
	}
*/	
	return score;
}




