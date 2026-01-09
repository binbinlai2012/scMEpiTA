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

void parse_bedname( string filename, string &sample, string &type )
{
    vector<string > ps = parse_string( filename, '/' );
    string ps0 = ps[ps.size()-1];
	vector<string > ps1 = parse_string( ps0, '.' );
	if ( (int)ps1.size() !=3 )
	{
		cout<<"unexpected bedfile name "<<filename<<" in parse_sample step1 "<<endl;
		exit(1);
	}
	sample = ps1[0];
	type = ps1[1];

}

void readinbed( string infile, string sample, string type,
	map<string, map<string, int > > &bc_sample_count,
	map<string, map<string, map<string, int > > > &bc_tag_count, 
    set<string> &tags )
{
	
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
//    string sample = "";
//	string type = "";
//	parse_bedname( infile, sample, type );
    string line = "";
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
		vector<string > parseditem = parse_string( line );
		if ( parseditem.size() != 6 )
		{
			cout<<"error bed line: "<<line<<endl;
			exit(1);
		}
	
		string chr = parseditem[0];
		
		int start = atoi(parseditem[1].c_str() );
		int end = atoi(parseditem[2].c_str() );

		string name = parseditem[3];
		string bc = parse_bc( name );
		if ( type == "atac" )
		{
			if ( bc_sample_count.find(bc) != bc_sample_count.end() )
			{
				if ( bc_sample_count[bc].find( sample ) != bc_sample_count[bc].end() )
				{
					bc_sample_count[bc][sample] += 1;
				} else
				{
					bc_sample_count[bc][sample] = 1;
				}
			} else
			{
				bc_sample_count[bc][sample] = 1;
			}
		} else
        {
            tags.insert( type );
            if ( bc_tag_count.find( bc ) != bc_tag_count.end() )
            {
                if ( bc_tag_count[bc].find( sample ) != bc_tag_count[bc].end() )
                {
                    if ( bc_tag_count[bc][sample].find( type ) != bc_tag_count[bc][sample].end() )
                        bc_tag_count[bc][sample][type] += 1;
                    else
                        bc_tag_count[bc][sample][type] = 1;
                } else
                    bc_tag_count[bc][sample][type] = 1;
            } else
                bc_tag_count[bc][sample][type] = 1;
        }

	}
    inf.close();
}

void readin_dnafiles( string infile, 
	map<string, map<string, int > > &bc_sample_count,
	map<string, map<string, map<string, int > > > &bc_tag_count,
    set<string > &tags )
{
    ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;

    while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line );
        if ( ps.size() != 3 )
        {
            cout<<"Unexpected line in bedfile: "<<line<<endl; exit(1);
        }
        string sample = ps[0];
        string tag = ps[1];
        string filename = ps[2];
        readinbed( filename, sample, tag, bc_sample_count, bc_tag_count, tags );
    }
    inf.close();

}

void filter_dominant_samples( map<string, map<string, int> > &bcs_samples,
    map<string, pair<int, double > > &bcs_collision_rate_sample,
    map<string, multimap<int, string, greater<int > > > &bcs_sorted_samples,
    set<string > &filtered_bcs_dna,
    int thr )
{
    for ( map<string, map<string, int> >::iterator ite = bcs_samples.begin(); ite != bcs_samples.end(); ++ite )
    {
        string bc = ite->first;
        multimap<int, string, greater<int > > sorted_samples;
        int total_bc_counts = 0;
        for ( map<string, int>::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            sorted_samples.insert( make_pair( si->second, si->first ) );
            total_bc_counts += si->second;
            if ( si->second >= thr )
            {
                filtered_bcs_dna.insert( bc );
            }
        }
        double rate = (sorted_samples.begin()->first*1.0)/total_bc_counts;
        bcs_collision_rate_sample.insert( make_pair( bc, make_pair(total_bc_counts, rate ) ) );
        bcs_sorted_samples.insert( make_pair(bc, sorted_samples ) );

    }
}

// cb_gc_uc.stat.txt generated by featureCount
void readinRNA1( string infile, map<string, int > &bc_rna_count )
{
    ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    string line;
    getline(inf, line);
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
        vector<string > parseditem = parse_string( line );
		if ( parseditem.size() != 4 )
		{
			cout<<"error bed line: "<<line<<endl;
			exit(1);
		}
        string bc = parseditem[1];
        int umic = atoi( parseditem[3].c_str() );
        bc_rna_count.insert( make_pair(bc, umic ) );

    }
    inf.close();
}

void readinRNA2( string infile, map<string, int > &bc_rna_count )
{
    ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
    string line;
    getline(inf, line);
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
        vector<string > parseditem = parse_string( line );
		if ( parseditem.size() != 5 )
		{
			cout<<"error bed line: "<<line<<endl;
			exit(1);
		}
        string bcname = parseditem[0];
        vector<string > ps2 = parse_string( bcname, ':' );
        string bc = ps2[ps2.size()-1];
        bc = bc.substr(0, bc.size()-1);
        int umic = atoi( parseditem[1].c_str() );
        bc_rna_count.insert( make_pair(bc, umic ) );

    }
    inf.close();
}

void filter_bcs_RNA( map<string, int > &bc_rna_count, set<string > &filtered_bcs, int thr )
{
    for ( map<string, int >::iterator ite = bc_rna_count.begin(); ite != bc_rna_count.end(); ++ite )
    {
        if ( ite->second > thr )
        {
            filtered_bcs.insert( ite->first );
        }
    }
}

void output( string outfile, 
    map<string, pair<int, double > > &bcs_collision_rate_sample,
    map<string, multimap<int, string, greater<int > > > &bcs_sorted_samples,
    set<string > &filtered_bcs_dna,
    map<string, map<string, map<string, int > > > &bc_tag_count,
    set<string > &tags,
    map<string, int > &bc_rna_count1,
    map<string, int > &bc_rna_count2,
    set<string > &filtered_bcs_rna,
    int thr_sample )
{

    set<string > filtered_bcs = filtered_bcs_dna;
    filtered_bcs.insert( filtered_bcs_rna.begin(), filtered_bcs_rna.end() );


    ofstream outf( outfile.data() );
    outf<<"barcodes\tUMIC1\tUMIC2\tATAC_count\tCollision_Rate\tSamples";
    for ( set<string >::iterator ite = tags.begin(); ite != tags.end(); ++ite )
    {
        outf<<"\t"<<*ite;
    }
    outf<<endl;
    for ( set<string>::iterator ite = filtered_bcs.begin(); ite != filtered_bcs.end(); ++ite )
    {
        string bc = *ite;
        outf<<bc;
        int umic1 = 0;
        if ( bc_rna_count1.find( bc ) !=  bc_rna_count1.end() )
            umic1 = bc_rna_count1[bc];
        outf<<"\t"<<umic1;

        int umic2 = 0;
        if ( bc_rna_count2.find( bc ) !=  bc_rna_count2.end() )
            umic2 = bc_rna_count2[bc];
        outf<<"\t"<<umic2;

        pair<int, double > sc = make_pair(0,0);
        if ( bcs_collision_rate_sample.find( bc ) != bcs_collision_rate_sample.end() )
            sc = bcs_collision_rate_sample[bc];
        outf<<"\t"<<sc.first<<"\t"<<sc.second;

        string dominant_sample = "";
        if ( bcs_sorted_samples.find( bc ) != bcs_sorted_samples.end() )
        {
            bool p = false;
            
            for ( multimap<int, string, greater<int > >::iterator si = bcs_sorted_samples[bc].begin(); 
                si != bcs_sorted_samples[bc].end(); ++si )
            {
                if ( si == bcs_sorted_samples[bc].begin() )
                {
                    outf<<"\t";
                    outf<<si->second<<":"<<si->first;
                    p = true;
                    dominant_sample = si->second;
                } else if ( si->first >= thr_sample )
                {
                    outf<<"|";
                    outf<<si->second<<":"<<si->first;
                    p = true;
                   
                } else
                    break;

            }
            if ( !p )
            {
                outf<<"\tNone:0";
            }
        } else
        {
            outf<<"\tNone:0";
        }

        for ( set<string >::iterator ti = tags.begin(); ti != tags.end(); ++ti )
        {
            string tg = *ti;
            bool fd = false;
            if ( bc_tag_count.find( bc ) != bc_tag_count.end() )
            {
                if ( bc_tag_count[bc].find( dominant_sample ) != bc_tag_count[bc].end() )
                {
                    if ( bc_tag_count[bc][dominant_sample].find( tg) !=  bc_tag_count[bc][dominant_sample].end() )
                    {
                        int c = bc_tag_count[bc][dominant_sample][tg];
                        fd = true;
                        outf<<"\t"<<c;
                    }

                }
            }
            if ( !fd )
                outf<<"\t0";

        }
        outf<<endl;

    }
    outf.close();

}

void exit_with_help()
{

	cerr <<"Usage:	Align_dna_rna_barcodes [OPTION1] [VALUE1] [[OPTION2] [VALUE2] ...]" <<endl;
	cerr <<"Options:" <<endl;
	
	cerr <<"-d   input files with each line a dna bed filename (Required)"<<endl;
	cerr <<"-r   RNA output file by featureCount (Required)" <<endl;
    cerr <<"-R   RNA output file by velocyto (Optional)"<<endl;
	cerr <<"-o   output file" <<endl;
    cerr <<"-t   filter threshold of dna count for single sample for barcodes (defalt 200)"<<endl;
	cerr <<"-T   filter threshold of umi count for barcodes (defalt 200)"<<endl;
	
	exit(1);
}
void exit_with_help(const char error[])
{
	cerr <<"Error:	" <<error <<endl;
	exit_with_help();

	exit(1);
}

int main( int argc, char* argv[] )
{
    string input_dnafiles = "";
    string input_rnafile1 = "";
    string input_rnafile2 = "";
    int thr_dna = 200;
    int thr_rna = 200;
    string outfile = "out.aligned_dna_rna_bc.txt";
    if (argc == 1)
	{
		exit_with_help();
	}
    for(int i=1; i<argc; i++)
	{
		if(argv[i][0] != '-')
			exit_with_help("Options must start with \'-\'.");

		if(argv[i][2] != '\0')
			exit_with_help("The option should be exactly one letter.");
		int option = argv[i][1];

		i++;

		switch(option)
		{
		
		case 'd':
			input_dnafiles = argv[i];
			break;
        case 'r':
            input_rnafile1 = argv[i];
            break;
        case 'R':
            input_rnafile2 = argv[i];
            break;
        case 't':
            thr_dna =  atoi( argv[i] );
            break;
        case 'T':
            thr_rna =  atoi( argv[i] );
            break;
        case 'o':
            outfile =  argv[i];
            break;
        default:
			exit_with_help();
		}
	}

    if ( input_dnafiles.empty() )
    {
        exit_with_help("Error: input_dnafiles empty!");
    }
    if ( input_rnafile1.empty() )
    {
        exit_with_help("Error: input_rnafile1 empty!");
    }

    map<string, map<string, int > > bc_sample_count;
	map<string, map<string, map<string, int > > > bc_tag_count;
    set<string> tags;
    cout<<"readin dna bed files"<<endl;
    readin_dnafiles( input_dnafiles, bc_sample_count, bc_tag_count, tags );

    map<string, pair<int, double > > bcs_collision_rate_sample;
    map<string, multimap<int, string, greater<int > > > bcs_sorted_samples;
    set<string > filtered_bcs_dna;
    cout<<"filter dominant dna bc"<<endl;
    filter_dominant_samples( bc_sample_count, bcs_collision_rate_sample, bcs_sorted_samples, filtered_bcs_dna, thr_dna );

    map<string, int > bc_rna_count1 ;
    cout<<"read in rna file 1"<<endl;
    readinRNA1( input_rnafile1, bc_rna_count1 );

    set<string > filtered_bcs_rna;
    filter_bcs_RNA( bc_rna_count1, filtered_bcs_rna, thr_rna );

    map<string, int > bc_rna_count2;
    if ( !input_rnafile2.empty() )
    {
        readinRNA2( input_rnafile2, bc_rna_count2 );
    }

 //   if ( outprefix[outprefix.length() - 1] != '/' )
//        outprefix+= '_';
 //   string outfile = outprefix+"aligned_dna_rna_bc.txt";
    output( outfile, 
        bcs_collision_rate_sample,
        bcs_sorted_samples,
        filtered_bcs_dna,
        bc_tag_count,
        tags,
        bc_rna_count1,
        bc_rna_count2,
        filtered_bcs_rna,
        thr_dna );

    return 1;
}

