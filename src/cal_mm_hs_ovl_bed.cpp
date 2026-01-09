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

void parse_bc( string bc, string& cell, string &header)
{
    vector<string > ps = parse_string( bc, '_' );
    if ( ps.size() != 3 )
    {
        cout<<"error unexpected bc "<<bc<<endl; exit(1);

    }
    cell = ps[1];
    header = ps[0];
}

void readinbed( string infile, map<string, set<string > > &cb_header )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;

    long k = 0;
    while(!inf.eof())
	{
        if ( k % 1000000 == 0 )
        {
            cout<<k;
        }
        k+= 1;

		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line );
        string bc = ps[3];
    //    bcs.insert( bc );
        string cell = "";
        string header = "";
        parse_bc( bc, cell, header );
        cb_header[cell].insert( header );
        
    }
    inf.close();

}


void cl_bed( string infile, map<string, set<string > > &cb_header )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;

	string outfile1 = infile+".ove.bed";
	string outfile2 = infile+".uni.bed";
	ofstream outf1(outfile1.data());
	ofstream outf2( outfile2.data() );
	string line;

    long k = 0;
    while(!inf.eof())
	{
        if ( k % 1000000 == 0 )
        {
            cout<<k;
        }
        k+= 1;
        
		getline(inf,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line );
        string bc = ps[3];
    //    bcs.insert( bc );
        string cell = "";
        string header = "";
        parse_bc( bc, cell, header );

        bool fd = false;
        if ( cb_header.find( cell ) != cb_header.end() )
        {
            if ( cb_header[cell].find( header ) != cb_header[cell].end() )
            {
                fd = true;
            }
        }
        if ( fd )
            outf1<<line<<endl;
        else
            outf2<<line<<endl;
    }
    inf.close();
	outf1.close();
	outf2.close();
}

int main( int argc, char* argv[] )
{
    if ( argc == 1 )
    {
	    cout<<"check overlap bcs from inbed_2 that align to inbed_1 "<<endl;
        cout<<"Usage: prog inbed_1 inbed_2"<<endl;
        exit(1);
    }
	string bedfile_1 = argv[1];
	string bedfile_2 = argv[2];
	
    map<string, set<string > > cb_header;
	readinbed( bedfile_1, cb_header );

	cl_bed( bedfile_2, cb_header );
	
	return 1;
}
