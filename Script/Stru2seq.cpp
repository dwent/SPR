//== Usage:  ./Stru2Seq  ../Domain.ID  ../SingleDomain/  ../Sequence/

#include <vector>
#include <fstream>
#include <iostream>
#include <string>


using namespace std;

/**********************************************************************
>> Translate residue names from 3-lette to 1-letter form
>> 
**********************************************************************/
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
		name1='?';
		cout << "\nAAname3to1()--Error: no matched residue name.\t" << name3 << endl;
	}
	
	return name1;
}

int main(int *argc, char *argv[])
{
	int i=0, j=0, k=0;
	
	vector<string> origID;
	vector<string> chainID;
	
	
	ifstream fin1(argv[1]);
	if( !fin1.is_open() )
	{
		cout << "\nError: can not open " << argv[1] << endl;
		exit(-1);
	}
	
	string stemp="";
	while( !fin1.eof() )
	{
		stemp="";
		fin1 >> stemp;
		if( stemp.size()>1 )
		{
			origID.push_back( stemp );
			//chainID.push_back( stemp.substr(4,1) );
		}
	}
	fin1.close();
	
	for(i=0; i<origID.size(); i++)
	{
		stemp="";
		stemp=origID[i];

		cout << i+1 << '\t' << stemp << endl;


		char chtemp[1000]="";
		strcpy(chtemp, argv[2]);
		strcat(chtemp, origID[i].c_str());
		strcat(chtemp, ".ent");
		
		
		ifstream fin(chtemp);
		if( !fin.is_open() )
		{
			cout << "\nError: can not open " << chtemp << endl;
			fin.close();
			continue;
		}
		
		
		char outname[1000]="";
		strcpy(outname, argv[3]);
		strcat(outname, stemp.c_str());
		strcat(outname, ".fasta");
		
		ofstream fout(outname);
		if( !fout.is_open() )
		{
			cout << "\nError: can not written into " << outname << endl;
			fout.close();
			continue;
		}
		fout << "> " << stemp << endl;
		
		
		//const char *c_chainid=chainID[i].c_str();
		char buf[100];
		bool MultiModel=false;
		
		string seq;
		
		while( !fin.eof() )
		{
			fin.getline(buf, 100);
			
			string sbuf(buf);
			
			if((sbuf.substr(0,5)=="MODEL")&&MultiModel) break;//MODEL 2 start
			else if(sbuf.substr(0,5)=="MODEL") MultiModel=true;
			
			if(strncmp(buf, "ATOM", 4)==0)
			{
				if(sbuf.substr(12,4)!=" CA ") continue;
			
				if((buf[16]!=' ')&&(buf[16]!='A')) continue;//retain A atom if the atom can not been determined
					
				string resid=sbuf.substr(17,3);
				
				seq.push_back(AAname3to1(resid));
			}
			
			
		}
		
		fout<<seq<<endl;
		
		fin.close();
		fout.close();
		
		
	}
	
	
	
	return 0;
}
