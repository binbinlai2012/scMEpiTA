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

string inttostr( int i )
{
    stringstream ss;
    ss << i;
    
    string st = ss.str();
    return st;
}

void readin_cb( string infile, set<string > &bcs )
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
    }
    inf.close();

}

void readinpeaks( string infile, map<string, set<pair<int, int > > > &peaks )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
    while(!inf.eof())
	{
        
		getline(inf,line);
		if ( line.empty() )
			break;
		vector<string > ps = parse_string( line );
		string chr = ps[0];
		int st = atoi(ps[1].c_str() );
		int ed = atoi(ps[2].c_str() );

		peaks[chr].insert(make_pair(st, ed) );
	}
	inf.close();
}

string parse_bc( string name )
{
	vector<string > ps = parse_string( name, '_' );
	if ( (int)ps.size() != 2 )
	{
		cout<<"unexpected str "<<name<<" in parsing barcode" <<endl;
		exit(1);
	}
	return ps[1];
}

string cal_overlap( string chr, int st, int ed, map<string, set<pair<int, int > > > &peaks  )
{
    string p = "None";
    bool fd = false;
    for ( set<pair<int, int > >::iterator ite = peaks[chr].begin(); ite != peaks[chr].end(); ++ite )
    {
        if ( ite->second < st )
            continue;
        if ( ite->first > ed )
            break;
        fd = true;
        p = chr+"-"+inttostr(ite->first)+"-"+inttostr(ite->second);
        break;

    }
    return p;
}

void assign_read_to_peaks( string infile, 
	set<string > &bcs, 
	map<string, set<pair<int, int > > > &peaks, 
	map<string, map<string, int > > &bc_peak_count,
    map<string, int > &unmapped )
{
	ifstream inf( infile.data() );
	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;

	while(!inf.eof())
	{
        
		getline(inf,line);
		if ( line.empty() )
			break;
		vector<string > ps = parse_string( line, '\t' );
        string bc = parse_bc( ps[3]);
        if ( bcs.find(bc) == bcs.end() )
            continue;
        string chr = ps[0];
        
        int st = atoi(ps[1].c_str() );
		int ed = atoi(ps[2].c_str() );

        string ovlp = cal_overlap( chr, st, ed, peaks );
        if ( ovlp != "None")
        {
            if ( bc_peak_count[bc].find( ovlp ) != bc_peak_count[bc].end() )
            {
                bc_peak_count[bc][ovlp] += 1;
            } else
            {
                bc_peak_count[bc][ovlp] = 1;
            }
        } else
        {
            if ( unmapped.find( bc) != unmapped.end() )
            {
                unmapped[bc] += 1;
            } else
                unmapped[bc] = 1;
        }
	}
    inf.close();

}
void assign_read_to_peaks_b( string infile, 
	set<string > &bcs, 
	map<string, set<pair<int, int > > > &peaks, 
	map<string, map<string, int > > &bc_peak_count,
    map<string, int > &unmapped )
{
    ifstream inf( infile.data() );
	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;

	while(!inf.eof())
	{
        
		getline(inf,line);
		if ( line.empty() )
			break;
        assign_read_to_peaks( line, bcs, peaks, bc_peak_count, unmapped );
    }
    inf.close();
}

void output( string outprefix, 
    set<string > &bcs, 
    map<string, set<pair<int, int > > > &peaks,
    map<string, map<string, int > > &bc_peak_count,
    map<string, int > &unmapped )
{
    string outfile1 = outprefix+".barcodes.tsv";
    ofstream outf1( outfile1.data() );
    string outfile2 = outprefix+".peaks.bed";
    ofstream outf2( outfile2.data() );
    string outfile3 = outprefix+".mtx";
    ofstream outf3( outfile3.data() );

    vector<string > bc_ve;
    for ( set<string>::iterator ite = bcs.begin(); ite != bcs.end(); ++ite )
    {

        outf1<<*ite<<endl;
        bc_ve.push_back(*ite);
    }
    outf1.close();

    vector<string > peak_ve;

    for ( map<string, set<pair<int, int > > >::iterator ite = peaks.begin(); ite != peaks.end(); ++ite )
    {
        string chr = ite->first;
        for ( set<pair<int, int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {

            string p = chr+"-"+inttostr(si->first)+"-"+inttostr(si->second);
            peak_ve.push_back(p);
            outf2<<p<<endl;
        }
    }
    outf2.close();

    outf3<<"%%MatrixMarket matrix coordinate integer general"<<endl;
    int bc_total = (int)bc_ve.size();
    int peak_total = (int)peak_ve.size();
    int c_total = 0;
    for ( map<string, map<string, int > >::iterator ite = bc_peak_count.begin(); ite != bc_peak_count.end(); ++ite )
    {
        for ( map<string, int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            c_total+=si->second;
        }
    }
    outf3<<peak_total<<"\t"<<bc_total<<"\t"<<c_total<<endl;


    for ( size_t i = 0; i < bc_ve.size(); ++i )
    {
        string bc = bc_ve[i];
        for ( size_t j = 0; j < peak_ve.size(); ++j )
        {
            string p = peak_ve[j];
            if ( bc_peak_count[bc].find( p ) != bc_peak_count[bc].end() )
            {
                outf3<<(j+1)<<"\t"<<(i+1)<<"\t"<<bc_peak_count[bc][p]<<endl;
            }
        }
    }
    outf3.close();

    string outfile4 = outprefix+".out_of_peaks.txt";
    ofstream outf4( outfile4.data() );
    for ( map<string, int >::iterator ite = unmapped.begin(); ite != unmapped.end(); ++ite )
    {
        outf4<<ite->first<<"\t"<<ite->second<<endl;
    }
    outf4.close();
}

int main( int argc, char * argv[] )
{
    if ( argc != 5 )
    {
        cout<<"generate mtx for epigenomic data"<<endl;
        cout<<"Usge: prog barcodefile peak.bed bedfiles[list] outprefix"<<endl;
        exit(1);
    }
    string bcfile = argv[1];
    string peakfile = argv[2];
    string bedfile = argv[3];
    string outprefix = argv[4];

    set<string > bcs;
    readin_cb( bcfile, bcs );


    map<string, set<pair<int, int > > > peaks;
    readinpeaks( peakfile, peaks );

    map<string, map<string, int > > bc_peak_count;
    map<string, int > unmapped;
    assign_read_to_peaks( bedfile, bcs, peaks, bc_peak_count, unmapped );

    cout<<"output"<<endl;
    output( outprefix,  bcs,  peaks, bc_peak_count, unmapped );

    return 1;
}
