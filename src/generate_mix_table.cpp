#include <iostream>
#include <string>
#include <map>
#include <list>
#include <vector>
#include <set>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <cstdlib>
#include <random>

using namespace std;

vector<string > parse_string( string & instr, char spl )
{
	vector<string > strve;
	string s = "";
	for ( size_t i = 0; i < instr.size(); ++i )
	{
		if ( instr[i] == spl )
		{
			if ( !s.empty() )
			{
				strve.push_back(s);
				s = "";
			}
		} else
		{
			s += instr[i];
		}
	}
	if ( !s.empty() )
	{
		strve.push_back(s);
	}
	return strve;
}

vector<string > parse_string( string & instr)
{
	vector<string > strve;
	string s = "";
	for ( size_t i = 0; i < instr.size(); ++i )
	{
		if ( instr[i] == '\t' || instr[i] == ' ' )
		{
			if ( !s.empty() )
			{
				strve.push_back(s);
				s = "";
			}
		} else
		{
			s += instr[i];
		}
	}
	if ( !s.empty() )
	{
		strve.push_back(s);
	}
	return strve;
}

void readinfile( string infile, set<string > &bcs, map<string, int > &cb_count_rna, map<string, int > &cb_count_atac )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
    getline(inf, line);
    while(!inf.eof())
	{
    
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line );
        string bc = ps[0];
        bcs.insert( bc );
        int cr = atoi( ps[1].c_str() );
        int ca = atoi( ps[3].c_str() );
        cb_count_rna.insert( make_pair(bc, cr ));
		cb_count_atac.insert( make_pair(bc, ca ) );
    }
    inf.close();


}

void output( string outfile, 
	set<string > &bcs, 
	map<string, int > &cb_count_rna_hg, 
	map<string, int > &cb_count_rna_mm, 
	map<string, int > &cb_count_atac_hg, 
	map<string, int > &cb_count_atac_mm )
{
	ofstream outf( outfile.data() );
	outf<<"barcode\tRNAHG\tRNAMM\tATACHG\tATACMM"<<endl;
	for ( set<string >::iterator ite = bcs.begin(); ite != bcs.end(); ++ite )
	{
		string bc = *ite;
		int rna_hg = 0;
		int rna_mm = 0;
		int atac_hg = 0;
		int atac_mm = 0;
		if ( cb_count_rna_hg.find( bc ) != cb_count_rna_hg.end() )
			rna_hg = cb_count_rna_hg[bc];
		if ( cb_count_rna_mm.find( bc ) != cb_count_rna_mm.end() )
			rna_mm = cb_count_rna_mm[bc];
		if ( cb_count_atac_hg.find( bc ) != cb_count_atac_hg.end() )
			atac_hg = cb_count_atac_hg[bc];
		if ( cb_count_atac_mm.find( bc ) != cb_count_atac_mm.end() )
			atac_mm = cb_count_atac_mm[bc];
		outf<<bc<<"\t"<<rna_hg<<"\t"<<rna_mm<<"\t"<<atac_hg<<"\t"<<atac_mm<<endl;
		
	}
	outf.close();
}

int main( int argc, char* argv[] )
{
	if ( argc == 1 )
	{
		cout<<"generate mix count table for HG and MM"<<endl;
		cout<<"Usage: prog HGinfile MMinfile outfile"<<endl;
		exit(1);
	}

	string infile_hg = argv[1];
	string infile_mm = argv[2];
	string outfile = argv[3];

	set<string > bcs;
	map<string, int > cb_count_rna_hg;
	map<string, int > cb_count_rna_mm;
	map<string, int > cb_count_atac_hg;
	map<string, int > cb_count_atac_mm;

	readinfile( infile_hg, bcs, cb_count_rna_hg, cb_count_atac_hg );
	readinfile( infile_mm, bcs, cb_count_rna_mm, cb_count_atac_mm );

	output( outfile, bcs, cb_count_rna_hg, cb_count_rna_mm, cb_count_atac_hg, cb_count_atac_mm );

	return 1;
}

