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

bool parse_line( string line, string &CB, string &UB,  string &gene )
{
    bool pass = false;
    vector<string > fds = parse_string( line );
    if ( fds.size() < 11 )
    {
        cout<<"unexpected sam line "<<fds.size()<<" "<<line<<endl; exit(1);
    }

    if ( fds.size() > 11 )
    {
        bool cb = false;
        bool ub = false;
        bool as = false;
        bool gn = false;
        for ( size_t i = 11; i < fds.size(); ++i )
        {
            vector<string > fd2 = parse_string( fds[i], ':' );
            if ( fd2.size() != 3 )
            {
                cout<<"unexpected tag field "<<fds[i]<<endl; exit(1);
            }
            string tag = fd2[0];
            if ( tag == "CB" )
            {
                CB = fd2[2];
                cb = true;
            } else if ( tag == "UB" )
            {
                ub = true;
                UB = fd2[2];
            } else if ( tag == "XS")
            {
                if ( fd2[2] == "Assigned")
                {
                    as = true;
                } 
            } else if ( tag == "XT" )
            {
                gene = fd2[2];
                gn = true;
            }
        }

        if ( cb && ub && as && gn )
        {
            pass = true;
        }
    }

    return pass;

}

void readinsam( string infile, map<string, map<string, set<string > > > &cb_gene_umi )
{
    ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;

    int i = 0;
	while(!inf.eof())
	{
        i+=1;
		getline(inf,line);
		if ( line.empty() )
			break;
        if ( line[0] == '@' )
            continue;
        if ( i %100000 == 0 )
            cout<<"processed "<<i<<" lines"<<endl;
        string CB = "";
        string UB = "";
        string gene = "";
        bool pass = parse_line( line, CB, UB,  gene );
        if ( pass )
        {
            cb_gene_umi[CB][gene].insert( UB );
        }
        
    }
    inf.close();

}

void cal( string outfile1,  map<string, map<string, set<string > > > &cb_gene_umi )
{
    ofstream outf1( outfile1.data() );

    map<string, pair<int, int > > cb_gc_uc;
    for ( map<string, map<string, set<string > > >::iterator ite = cb_gene_umi.begin(); ite != cb_gene_umi.end(); ++ite )
    {
        int gc = 0;
        int uc = 0;
        for ( map<string, set<string > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            gc += 1;
            uc += (int)si->second.size();
        }
        cb_gc_uc.insert( make_pair(ite->first, make_pair(gc, uc ) ) );
    }

    outf1<<"Rank\tCell_Barcode\tGeneCount\tUMICount"<<endl;
    multimap<int, string, greater<int> > sorted_cb;
    for ( map<string, pair<int, int > >::iterator ite = cb_gc_uc.begin(); ite != cb_gc_uc.end(); ++ite )
    {
        string cb = ite->first;
        int gc = ite->second.first;
        int uc = ite->second.second;
        sorted_cb.insert(make_pair(uc, cb ) );
    }
    int rk = 1;
    for ( multimap<int, string, greater<int> >::iterator ite = sorted_cb.begin(); ite != sorted_cb.end(); ++ite )
    {
        string cb = ite->second;
        int gc = cb_gc_uc[cb].first;
        int uc = cb_gc_uc[cb].second;
        outf1<<rk<<"\t"<<cb<<"\t"<<gc<<"\t"<<uc<<endl;
        rk+=1;
    }
    outf1.close();
}

int main(int argc, char* argv[] )
{
    if( argc == 1 )
    {
        cout<<"statistics of cellbarcode genecount umicount table "<<endl;
        cout<<"Usage: prog insam outfile"<<endl;
        exit(1);
    }

    string insam = argv[1];
    string outfile = argv[2];

    map<string, map<string, set<string > > > cb_gene_umi;
    readinsam( insam, cb_gene_umi );

    cal( outfile, cb_gene_umi );

    return 1;

}
