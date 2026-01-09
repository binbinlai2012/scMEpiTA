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

void func( string infile, string outfile1, string outfile2 )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;

	ofstream outf1( outfile1.data() );
	ofstream outf2( outfile2.data() );

    int i = 0;
	getline(inf, line);
	outf1<<line<<endl;
	outf2<<line<<endl;
	while(!inf.eof())
	{
        i+=1;
		getline(inf,line);
		if ( line.empty() )
			break;
        
		vector<string > ps = parse_string(line, '\t');
		if ( ps[0] != "@SQ")
			continue;
		
		vector<string > spchr1 = parse_string( ps[1], ':' );
		vector<string > spchr2 = parse_string( spchr1[1], '_');
		string sp = spchr2[0];
		string chr = spchr2[1];
		if ( sp == "GRCh38" )
		{
			for ( size_t i = 0; i < ps.size(); ++i )
			{
				if ( i > 0 )
					outf1<<"\t";
				if ( i == 1 )
					outf1<<"SN:"<<chr;
				else
					outf1<<ps[i];
			}
			outf1<<endl;
		} else if ( sp == "mm10" )
		{
			for ( size_t i = 0; i < ps.size(); ++i )
			{
				if ( i > 0 )
					outf2<<"\t";
				if ( i == 1 )
					outf2<<"SN:"<<chr;
				else
					outf2<<ps[i];
			}
			outf2<<endl;
		} 
	}
    outf1.close();
    outf2.close();
    inf.close();

}

int main( int argc, char* argv[] )
{
    if ( argc == 1 )
	{
		cout<<"split combined sam header"<<endl;
		cout<<"Usage: prog combsamheader outhgsamheader outmmsamheader"<<endl; exit(1);
	}

	string infile = argv[1];
	string outfile1 = argv[2];
	string outfile2 = argv[3];

	func( infile, outfile1, outfile2 );

	return 1;
}
