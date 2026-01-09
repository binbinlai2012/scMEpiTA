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

/*
string parse_sample( string name )
{
	vector<string > ps1 = parse_string( name, '.' );
	if ( (int)ps1.size() < 2 )
	{
		cout<<"unexpected file name "<<name<<" in parse_sample step1 "<<endl;
		exit(1);
	}
	vector<string > ps2 = parse_string( ps1[0], '-');
	if ( (int)ps2.size() < 2 )
	{
		cout<<"unexpected file name "<<name<<" in parse_sample step2 "<<endl;
		exit(1);
	}
	string sp = ps2[1];
	if ( (int)ps2.size()>2 )
	{
		for ( size_t i = 2; i < ps2.size(); ++i )
		{
			sp+="-";
			sp+=ps2[i];
		}
	}
	return sp;
} */

// For atac_R2 (both sc and sa libs), it is known which samples the reads (whole file) are from, which can be recognized by file name sa/sc*-$ATAC{j}.*.R2.bam (here input from bedfiles)
// For sc/sa libs, if the R2 are samples_index, the sample is also known for R1 (no matter R1 is atac or cuttag), which can be recognized by file name with patterns sa*-$ATAC{j}.*.R1.bam
// First proceed atac_R2. And construct sample_barcodes maps.
// Here barcodes may go to different samples, which indicates collision.
void read_bed_R2_sampleknown_atac( string &infile, string insample,
	map<string, map<string, int> > &bcs_samples,
	map<string, map<string, map< pair<int, int >, string > > > &sample_chr_pos_line_atac )
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
 //   string insample = parse_sample(infile);    // 
	string line;
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
	//	if ( line[0] == '#' )
	//		continue;
		
		vector<string > parseditem = parse_string( line );
		if ( parseditem.size() != 6 )
		{
			cout<<"error bed line: "<<line<<endl;
			exit(1);
		}
	
		string chr = parseditem[0];
		
		int start = atoi(parseditem[1].c_str() );
		int end = atoi(parseditem[2].c_str() );
		sample_chr_pos_line_atac[insample][chr][make_pair(start, end )] = line;

		string name = parseditem[3];
		string bc = parse_bc( name );

        if ( bcs_samples.find(bc) == bcs_samples.end() )
    		bcs_samples[bc].insert( make_pair(insample, 1 ) );
        else 
        {
            if ( bcs_samples[bc].find( insample ) != bcs_samples[bc].end() )
                bcs_samples[bc][insample] += 1;
            else
                bcs_samples[bc].insert( make_pair(insample, 1) );
        }
	}
	
	inf.close();
	
}

// for bcs that have multiple samples. we select dominant sample for inferring at later stages. 
// Mark bcs that do not have dominant samples.
void set_dominant_samples( map<string, map<string, int> > &bcs_samples,
    map<string, pair<int, double > > &bcs_collision_rate_sample,
    map<string, string > &bcs_dominant_samples,
    set<string > &collision_bcs_sample,
    double collision_thr )
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
        }
        double rate = (sorted_samples.begin()->first*1.0)/total_bc_counts;
        bcs_collision_rate_sample.insert( make_pair( bc, make_pair(total_bc_counts, rate ) ) );
        if ( rate >= collision_thr )
        {
            bcs_dominant_samples.insert( make_pair(bc, sorted_samples.begin()->second ) );
        } else
        {
            collision_bcs_sample.insert( bc );
        }

    }
}

// For CUTTag, we need to know not only sample attr but also the Tag type. For cuttag_R2, it can be inferred from file name.
// So we proceed cuttag_R2 files. Here we can construct bc tag maps.
void read_bed_R2_tagknown_cuttag( string infile, string intag,
    map<string, string > &bcs_dominant_samples,
    set<string > &collision_bcs_sample,
    map< string, map<string,int> > &bcs_tags,
    set<string > &unknown_bc,
    map<string, map<string, map<string, map< pair<int, int >, string > > > > &sample_tag_chr_pos_line_cuttag )
{
    ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
 //   string intag = parse_sample(infile);    // use the same func as sample, becase they are at the same pos.
	string line;
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
        string name = parseditem[3];
		string bc = parse_bc( name );

        string chr = parseditem[0];
		
		int start = atoi(parseditem[1].c_str() );
		int end = atoi(parseditem[2].c_str() );

        if ( bcs_dominant_samples.find(bc) != bcs_dominant_samples.end() )
		{
			
			string sample = bcs_dominant_samples[bc];
			sample_tag_chr_pos_line_cuttag[sample][intag][chr][make_pair(start, end )] = line;

            if ( bcs_tags.find( bc ) == bcs_tags.end() )
                bcs_tags[bc].insert( make_pair( intag, 1) );
            else
            {
                if ( bcs_tags[bc].find( intag ) != bcs_tags[bc].end() )
                    bcs_tags[bc][intag] += 1;
                else
                    bcs_tags[bc][intag] = 1;
            }
		} else if ( collision_bcs_sample.find( bc ) == collision_bcs_sample.end() )
        {   
            unknown_bc.insert( bc );
        }

    }
    inf.close();
}

// for bcs that have multiple tags. we select dominant tag for inferring at later stages. 
// Mark bcs that do not have dominant tags.
void set_dominant_tags( map<string, map<string, int> > &bcs_tags,
    map<string, pair<int, double > > &bcs_collision_rate_tag,
    map<string, string > &bcs_dominant_tags,
    set<string > &collision_bcs_tag,
    double collision_thr )
{

    for ( map<string, map<string, int> >::iterator ite = bcs_tags.begin(); ite != bcs_tags.end(); ++ite )
    {
        string bc = ite->first;
        multimap<int, string, greater<int > > sorted_tags;
        int total_bc_counts = 0;
        for ( map<string, int>::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            sorted_tags.insert( make_pair( si->second, si->first ) );
            total_bc_counts += si->second;
        }
        double rate = (sorted_tags.begin()->first*1.0)/total_bc_counts;
        bcs_collision_rate_tag.insert( make_pair( bc, make_pair(total_bc_counts, rate ) ) );
        if ( rate >= collision_thr )
        {
            bcs_dominant_tags.insert( make_pair(bc, sorted_tags.begin()->second ) );
        } else
        {
            collision_bcs_tag.insert( bc );
        }

    }
}

// Proceed atac R1 (sa) for atac_R2. Although it is can be inferred from bc_sample map. But it can also be inferred from file name.
void read_bed_R1_sampleknown_atac( string &infile, string insample,
    map<string, map<string, map< pair<int, int >, string > > > &sample_chr_pos_line_atac )
{
    ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
 //   string insample = parse_sample(infile);
	string line;
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

        
		sample_chr_pos_line_atac[insample][chr][make_pair(start, end )] = line;


    }
    inf.close();
}

// Proceed atac R1 (sa) for cuttag_R2. Infer sample from bc sample maps.
void read_bed_R1_sampleunknown_atac( string &infile, 
	map<string, string > &bcs_dominant_samples,
    set<string > &collision_bcs_sample,
	set<string > &unknown_bc,
	map<string, map<string, map< pair<int, int >, string > > > &sample_chr_pos_line_atac )
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
		if ( bcs_dominant_samples.find(bc) != bcs_dominant_samples.end() )
		{
			string sample = bcs_dominant_samples[bc];
			sample_chr_pos_line_atac[sample][chr][make_pair(start, end )] = line;

		} else
        {
            if ( collision_bcs_sample.find( bc ) == collision_bcs_sample.end() )
                unknown_bc.insert( bc );
        }
	}
    inf.close();
} 

// Proceed R1 cuttag (sc lib) for cuttag_R2. the tag type can be inferred from file name. the sample can be inferred from bc sample maps.
void read_bed_R1_tagknown_sampleunknown_cuttag( string infile, string intag,
    map<string, string > &bcs_dominant_samples,
    set<string > &collision_bcs_sample,
    set<string > &unknown_bc,
    map<string, map<string, map<string, map< pair<int, int >, string > > > > &sample_tag_chr_pos_line_cuttag )
{
    ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
 //   string intag = parse_sample(infile);    // use the same func as sample, becase they are at the same pos.
	string line;
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
        string name = parseditem[3];
		string bc = parse_bc( name );

        string chr = parseditem[0];
		
		int start = atoi(parseditem[1].c_str() );
		int end = atoi(parseditem[2].c_str() );
        if ( bcs_dominant_samples.find(bc) != bcs_dominant_samples.end() )
		{
			
			string sample = bcs_dominant_samples[bc];
			sample_tag_chr_pos_line_cuttag[sample][intag][chr][make_pair(start, end )] = line;
		} else
        {
            if ( collision_bcs_sample.find(bc) == collision_bcs_sample.end() )
                unknown_bc.insert( bc ); 
        }
    }
}

// Last, For sc lib and atac_R2, sample can be inferred from file name, tag type can be inferred from bc tag maps. 
void read_bed_R1_sampleknown_tagunkown_cuttag( string &infile, string insample,
    map<string, string > &bcs_dominant_tags,
    set<string > &collision_bcs_tag,
	set<string > &unknown_tag,
	map<string, map<string, map<string, map< pair<int, int >, string > > > > &sample_tag_chr_pos_line_cuttag )
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
 //   string insample = parse_sample(infile);    
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
		if ( bcs_dominant_tags.find(bc) != bcs_dominant_tags.end() )
		{
			string tag = bcs_dominant_tags[bc];
			sample_tag_chr_pos_line_cuttag[insample][tag][chr][make_pair(start, end )] = line;

		} else
        {
            if ( collision_bcs_tag.find(bc) == collision_bcs_tag.end() )
                unknown_tag.insert(bc);
        }
	}
    inf.close();
}

void readinfiles_main( string infile, 
    map<string, map<string, map< pair<int, int >, string > > > &sample_chr_pos_line_atac,
    map<string, map<string, map<string, map< pair<int, int >, string > > > > &sample_tag_chr_pos_line_cuttag,
    map<string, map<string, int > > &bc_sample_count,
    map<string, pair<int, double > > &bcs_collision_rate_sample,
    map<string, string > &bcs_dominant_samples,
    set<string > &collision_bcs_sample,
    map<string, map<string, int > > &bc_tag_count,
    map<string, pair<int, double > > &bcs_collision_rate_tag,
    map<string, string > &bcs_dominant_tags,
    set<string> &collision_bcs_tag,
    map<string, set<string > > &unknown_sample_map,
    map<string, set<string > > &unknown_tag_map,
    double collision_thr )
{
    ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
    map<string, map<string, string > > type_sampletag_filename;
    while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
		vector<string > parseditem = parse_string( line );
		if ( parseditem.size() != 3 )
		{
			cout<<"error metadata line: "<<line<<endl;
			exit(1);
		}
        string type = parseditem[0];
        type_sampletag_filename[type][parseditem[1]] = parseditem[2];
        
    }
    inf.close();

    // first sa_sample_R2.
    cout<<">Readin sa_sample_R2"<<endl;
    if ( type_sampletag_filename.find("sa_sample_R2") == type_sampletag_filename.end() )
    {
        cout<<"error cannot get sa_sample_R2 files"<<endl; exit(1);
    }
    for ( map<string, string >::iterator ite = type_sampletag_filename["sa_sample_R2"].begin(); 
        ite != type_sampletag_filename["sa_sample_R2"].end(); ++ite )
    {
        string sample = ite->first;
        string filename = ite->second;
        read_bed_R2_sampleknown_atac( filename, sample, bc_sample_count, sample_chr_pos_line_atac );
    }
    // second sc_sample_R2.
    cout<<">Readin sc_sample_R2"<<endl;
    if ( type_sampletag_filename.find("sc_sample_R2") == type_sampletag_filename.end() )
    {
        cout<<"error cannot get sc_sample_R2 files"<<endl; exit(1);
    }
    for ( map<string, string >::iterator ite = type_sampletag_filename["sc_sample_R2"].begin(); 
        ite != type_sampletag_filename["sc_sample_R2"].end(); ++ite )
    {
        string sample = ite->first;
        string filename = ite->second;
        read_bed_R2_sampleknown_atac( filename, sample, bc_sample_count, sample_chr_pos_line_atac );
    }
    // calculate collision_sample
    cout<<"set dominant samples"<<endl;
    set_dominant_samples( bc_sample_count, bcs_collision_rate_sample, bcs_dominant_samples, collision_bcs_sample, collision_thr );
    
    // third sc_cuttag_R2.
    cout<<">Readin sc_tag_R2"<<endl;
    if ( type_sampletag_filename.find("sc_tag_R2") == type_sampletag_filename.end() )
    {
        cout<<"error cannot get sc_tag_R2 files"<<endl; exit(1);
    }
    for ( map<string, string >::iterator ite = type_sampletag_filename["sc_tag_R2"].begin(); 
        ite != type_sampletag_filename["sc_tag_R2"].end(); ++ite )
    {
        string sample = ite->first;
        string filename = ite->second;
        set<string > unkown_sample;
        read_bed_R2_tagknown_cuttag( filename, sample,
            bcs_dominant_samples, 
            collision_bcs_sample, 
            bc_tag_count, 
            unkown_sample, 
            sample_tag_chr_pos_line_cuttag );
        unknown_sample_map.insert( make_pair(filename, unkown_sample) );
    }
    // fourth sa_cuttag_R2.
    cout<<">Readin sa_tag_R2"<<endl;
    if ( type_sampletag_filename.find("sa_tag_R2") == type_sampletag_filename.end() )
    {
        cout<<"error cannot get sa_tag_R2 files"<<endl; exit(1);
    }
    for ( map<string, string >::iterator ite = type_sampletag_filename["sa_tag_R2"].begin(); 
        ite != type_sampletag_filename["sa_tag_R2"].end(); ++ite )
    {
        string sample = ite->first;
        string filename = ite->second;
        set<string > unkown_sample;
        read_bed_R2_tagknown_cuttag( filename, sample, 
            bcs_dominant_samples, 
            collision_bcs_sample, 
            bc_tag_count, 
            unkown_sample, 
            sample_tag_chr_pos_line_cuttag );
        unknown_sample_map.insert( make_pair(filename, unkown_sample) );
    }
    // calculate collision_tag
    cout<<"set dominant tags"<<endl;
    set_dominant_tags( bc_tag_count, bcs_collision_rate_tag, bcs_dominant_tags, collision_bcs_tag, collision_thr );
    
    // fifth: sa_sample_R1
    cout<<">Readin sa_sample_R1"<<endl;
    if ( type_sampletag_filename.find("sa_sample_R1") == type_sampletag_filename.end() )
    {
        cout<<"error cannot get sa_sample_R1 files"<<endl; exit(1);
    }
    for ( map<string, string >::iterator ite = type_sampletag_filename["sa_sample_R1"].begin(); 
        ite != type_sampletag_filename["sa_sample_R1"].end(); ++ite )
    {
        string sample = ite->first;
        string filename = ite->second;
        read_bed_R1_sampleknown_atac(filename, sample, sample_chr_pos_line_atac );
    }
    // sixth, sa_tag_R1
    cout<<">Readin sa_tag_R1"<<endl;
    if ( type_sampletag_filename.find("sa_tag_R1") == type_sampletag_filename.end() )
    {
        cout<<"error cannot get sa_tag_R1 files"<<endl; exit(1);
    }
    for ( map<string, string >::iterator ite = type_sampletag_filename["sa_tag_R1"].begin(); 
        ite != type_sampletag_filename["sa_tag_R1"].end(); ++ite )
    {
        string sample = ite->first;
        string filename = ite->second;
        set<string > unkown_sample;
        read_bed_R1_sampleunknown_atac( filename, 
            bcs_dominant_samples, 
            collision_bcs_sample, 
            unkown_sample, 
            sample_chr_pos_line_atac );
        unknown_sample_map.insert( make_pair(filename, unkown_sample) );
    }
    // 7th, sc_tag_R1
    cout<<">Readin sc_tag_R1"<<endl;
    if ( type_sampletag_filename.find("sc_tag_R1") == type_sampletag_filename.end() )
    {
        cout<<"error cannot get sc_tag_R1 files"<<endl; exit(1);
    }
    for ( map<string, string >::iterator ite = type_sampletag_filename["sc_tag_R1"].begin(); 
        ite != type_sampletag_filename["sc_tag_R1"].end(); ++ite )
    {
        string tag = ite->first;
        string filename = ite->second;
        set<string > unknown_sample;
        read_bed_R1_tagknown_sampleunknown_cuttag( filename, tag,
            bcs_dominant_samples, 
            collision_bcs_sample, 
            unknown_sample, 
            sample_tag_chr_pos_line_cuttag );
        unknown_sample_map.insert( make_pair(filename, unknown_sample) );
    }
    // 8th, sc_sample_R1
    cout<<">Readin sc_sample_R1"<<endl;
    if ( type_sampletag_filename.find("sc_sample_R1") == type_sampletag_filename.end() )
    {
        cout<<"error cannot get sc_sample_R1 files"<<endl; exit(1);
    }
    for ( map<string, string >::iterator ite = type_sampletag_filename["sc_sample_R1"].begin(); 
        ite != type_sampletag_filename["sc_sample_R1"].end(); ++ite )
    {
        string sample = ite->first;
        string filename = ite->second;
        set<string> unknown_tag;
        read_bed_R1_sampleknown_tagunkown_cuttag( filename, sample,
            bcs_dominant_tags, 
            collision_bcs_tag, 
            unknown_tag, 
            sample_tag_chr_pos_line_cuttag );
        unknown_tag_map.insert( make_pair(filename, unknown_tag ) );

    }


}

void readinfiles_main_short( string infile, 
    map<string, map<string, map< pair<int, int >, string > > > &sample_chr_pos_line_atac,
    map<string, map<string, map<string, map< pair<int, int >, string > > > > &sample_tag_chr_pos_line_cuttag,
    map<string, map<string, int > > &bc_sample_count,
    map<string, pair<int, double > > &bcs_collision_rate_sample,
    map<string, string > &bcs_dominant_samples,
    set<string > &collision_bcs_sample,
    map<string, map<string, int > > &bc_tag_count,
    map<string, pair<int, double > > &bcs_collision_rate_tag,
    map<string, string > &bcs_dominant_tags,
    set<string> &collision_bcs_tag,
    map<string, set<string > > &unknown_sample_map,
    map<string, set<string > > &unknown_tag_map,
    double collision_thr )
{
    ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
    map<string, map<string, string > > type_sampletag_filename;
    while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
		vector<string > parseditem = parse_string( line );
		if ( parseditem.size() != 3 )
		{
			cout<<"error metadata line: "<<line<<endl;
			exit(1);
		}
        string type = parseditem[0];
        type_sampletag_filename[type][parseditem[1]] = parseditem[2];
        
    }
    inf.close();

    // first sa_sample_R2.
    cout<<">Readin sa_sample_R2"<<endl;
    if ( type_sampletag_filename.find("sa_sample_R2") == type_sampletag_filename.end() )
    {
        cout<<"error cannot get sa_sample_R2 files"<<endl; exit(1);
    }
    for ( map<string, string >::iterator ite = type_sampletag_filename["sa_sample_R2"].begin(); 
        ite != type_sampletag_filename["sa_sample_R2"].end(); ++ite )
    {
        string sample = ite->first;
        string filename = ite->second;
        read_bed_R2_sampleknown_atac( filename, sample, bc_sample_count, sample_chr_pos_line_atac );
    }
    // second sc_sample_R2.
/*    cout<<">Readin sc_sample_R2"<<endl;
    if ( type_sampletag_filename.find("sc_sample_R2") == type_sampletag_filename.end() )
    {
        cout<<"error cannot get sc_sample_R2 files"<<endl; exit(1);
    }
    for ( map<string, string >::iterator ite = type_sampletag_filename["sc_sample_R2"].begin(); 
        ite != type_sampletag_filename["sc_sample_R2"].end(); ++ite )
    {
        string sample = ite->first;
        string filename = ite->second;
        read_bed_R2_sampleknown_atac( filename, sample, bc_sample_count, sample_chr_pos_line_atac );
    } */
    // calculate collision_sample
    cout<<"set dominant samples"<<endl;
    set_dominant_samples( bc_sample_count, bcs_collision_rate_sample, bcs_dominant_samples, collision_bcs_sample, collision_thr );
    
    // third sc_cuttag_R2.
/*    cout<<">Readin sc_tag_R2"<<endl;
    if ( type_sampletag_filename.find("sc_tag_R2") == type_sampletag_filename.end() )
    {
        cout<<"error cannot get sc_tag_R2 files"<<endl; exit(1);
    }
    for ( map<string, string >::iterator ite = type_sampletag_filename["sc_tag_R2"].begin(); 
        ite != type_sampletag_filename["sc_tag_R2"].end(); ++ite )
    {
        string sample = ite->first;
        string filename = ite->second;
        set<string > unkown_sample;
        read_bed_R2_tagknown_cuttag( filename, sample,
            bcs_dominant_samples, 
            collision_bcs_sample, 
            bc_tag_count, 
            unkown_sample, 
            sample_tag_chr_pos_line_cuttag );
        unknown_sample_map.insert( make_pair(filename, unkown_sample) );
    } */
    // fourth sa_cuttag_R2.
    cout<<">Readin sa_tag_R2"<<endl;
    if ( type_sampletag_filename.find("sa_tag_R2") == type_sampletag_filename.end() )
    {
        cout<<"error cannot get sa_tag_R2 files"<<endl; exit(1);
    }
    for ( map<string, string >::iterator ite = type_sampletag_filename["sa_tag_R2"].begin(); 
        ite != type_sampletag_filename["sa_tag_R2"].end(); ++ite )
    {
        string sample = ite->first;
        string filename = ite->second;
        set<string > unkown_sample;
        read_bed_R2_tagknown_cuttag( filename, sample, 
            bcs_dominant_samples, 
            collision_bcs_sample, 
            bc_tag_count, 
            unkown_sample, 
            sample_tag_chr_pos_line_cuttag );
        unknown_sample_map.insert( make_pair(filename, unkown_sample) );
    }
    // calculate collision_tag
    cout<<"set dominant tags"<<endl;
    set_dominant_tags( bc_tag_count, bcs_collision_rate_tag, bcs_dominant_tags, collision_bcs_tag, collision_thr );
    
    // fifth: sa_sample_R1
 /*   cout<<">Readin sa_sample_R1"<<endl;
    if ( type_sampletag_filename.find("sa_sample_R1") == type_sampletag_filename.end() )
    {
        cout<<"error cannot get sa_sample_R1 files"<<endl; exit(1);
    }
    for ( map<string, string >::iterator ite = type_sampletag_filename["sa_sample_R1"].begin(); 
        ite != type_sampletag_filename["sa_sample_R1"].end(); ++ite )
    {
        string sample = ite->first;
        string filename = ite->second;
        read_bed_R1_sampleknown_atac(filename, sample, sample_chr_pos_line_atac );
    }
    // sixth, sa_tag_R1
    cout<<">Readin sa_tag_R1"<<endl;
    if ( type_sampletag_filename.find("sa_tag_R1") == type_sampletag_filename.end() )
    {
        cout<<"error cannot get sa_tag_R1 files"<<endl; exit(1);
    }
    for ( map<string, string >::iterator ite = type_sampletag_filename["sa_tag_R1"].begin(); 
        ite != type_sampletag_filename["sa_tag_R1"].end(); ++ite )
    {
        string sample = ite->first;
        string filename = ite->second;
        set<string > unkown_sample;
        read_bed_R1_sampleunknown_atac( filename, 
            bcs_dominant_samples, 
            collision_bcs_sample, 
            unkown_sample, 
            sample_chr_pos_line_atac );
        unknown_sample_map.insert( make_pair(filename, unkown_sample) );
    }
    // 7th, sc_tag_R1
    cout<<">Readin sc_tag_R1"<<endl;
    if ( type_sampletag_filename.find("sc_tag_R1") == type_sampletag_filename.end() )
    {
        cout<<"error cannot get sc_tag_R1 files"<<endl; exit(1);
    }
    for ( map<string, string >::iterator ite = type_sampletag_filename["sc_tag_R1"].begin(); 
        ite != type_sampletag_filename["sc_tag_R1"].end(); ++ite )
    {
        string tag = ite->first;
        string filename = ite->second;
        set<string > unknown_sample;
        read_bed_R1_tagknown_sampleunknown_cuttag( filename, tag,
            bcs_dominant_samples, 
            collision_bcs_sample, 
            unknown_sample, 
            sample_tag_chr_pos_line_cuttag );
        unknown_sample_map.insert( make_pair(filename, unknown_sample) );
    }
    // 8th, sc_sample_R1
    cout<<">Readin sc_sample_R1"<<endl;
    if ( type_sampletag_filename.find("sc_sample_R1") == type_sampletag_filename.end() )
    {
        cout<<"error cannot get sc_sample_R1 files"<<endl; exit(1);
    }
    for ( map<string, string >::iterator ite = type_sampletag_filename["sc_sample_R1"].begin(); 
        ite != type_sampletag_filename["sc_sample_R1"].end(); ++ite )
    {
        string sample = ite->first;
        string filename = ite->second;
        set<string> unknown_tag;
        read_bed_R1_sampleknown_tagunkown_cuttag( filename, sample,
            bcs_dominant_tags, 
            collision_bcs_tag, 
            unknown_tag, 
            sample_tag_chr_pos_line_cuttag );
        unknown_tag_map.insert( make_pair(filename, unknown_tag ) );

    }

*/
}

void outputbed( string outfile, map<string, map< pair<int, int >, string > > &chr_pos_line )
{
    ofstream outf(outfile.data());
    for ( map<string, map< pair<int, int >, string > >::iterator ite = chr_pos_line.begin(); ite != chr_pos_line.end(); ++ite )
    {

        for ( map< pair<int, int >, string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            outf<<si->second<<endl;
        }
    }
    outf.close();
}

void output( map<string, map<string, map< pair<int, int >, string > > > &sample_chr_pos_line_atac,
    map<string, map<string, map<string, map< pair<int, int >, string > > > > &sample_tag_chr_pos_line_cuttag,
    map<string, map<string, int > > &bc_sample_count,
    map<string, pair<int, double > > &bcs_collision_rate_sample,
    map<string, string > &bcs_dominant_samples,
    set<string > &collision_bcs_sample,
    map<string, map<string, int > > &bc_tag_count,
    map<string, pair<int, double > > &bcs_collision_rate_tag,
    map<string, string > &bcs_dominant_tags,
    set<string> &collision_bcs_tag,
    map<string, set<string > > &unknown_sample_map,
    map<string, set<string > > &unknown_tag_map,
    string outprefix )
{
    if ( outprefix[outprefix.length() - 1] != '/' )
        outprefix+= '_';

    // output atac
    cout<<"output atac"<<endl;
    for ( map<string, map<string, map< pair<int, int >, string > > >::iterator ite = sample_chr_pos_line_atac.begin();
        ite != sample_chr_pos_line_atac.end(); ++ite )
    {
        string sample = ite->first;
        string outfile = outprefix + sample+".atac.bed";
        
        outputbed( outfile, ite->second );
    }

    // output tag
    cout<<"output tag"<<endl;
    for ( map<string, map<string, map<string, map< pair<int, int >, string > > > >::iterator ite = sample_tag_chr_pos_line_cuttag.begin();
        ite != sample_tag_chr_pos_line_cuttag.end(); ++ite )
    {
        string sample = ite->first;
        for ( map<string, map<string, map< pair<int, int >, string > > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string tag = si->first;
            string outfile = outprefix+sample+"."+tag+".bed";
            outputbed( outfile, si->second );

        }
    }

    // output barcode_sample_mapping
    cout<<"output barcode_sample_mapping"<<endl;
    string outfile_bc_sample = outprefix+"sample_barcodes_mapping.txt";
    ofstream outf_bs( outfile_bc_sample.data() );
    outf_bs<< "barcodes\tcounts\tcollision_rate\tdominant_sample"<<endl;
    for ( map<string, pair<int, double > >::iterator ite = bcs_collision_rate_sample.begin(); ite != bcs_collision_rate_sample.end(); ++ite )
    {
        string bc = ite->first;
        int count = ite->second.first;
        double rate = ite->second.second;
        
        if ( bcs_dominant_samples.find(bc) != bcs_dominant_samples.end() )
        {
            string ds = bcs_dominant_samples[bc];
            outf_bs<<bc<<"\t"<<count<<"\t"<<rate<<"\t"<<ds<<endl;
        } else
            outf_bs<<bc<<"\t"<<count<<"\t"<<rate<<"\tCollision"<<endl;
        
    }
    outf_bs.close();

    // output barcode_tag_mapping
    cout<<"output barcode_tag_mapping"<<endl;
    string outfile_bc_tag = outprefix+"tag_barcodes_mapping.txt";
    ofstream outf_bt( outfile_bc_tag.data() );
    outf_bt<<"barcodes\tcounts\tcollision_rate\tdominant_tag"<<endl;
    for ( map<string, pair<int, double > >::iterator ite = bcs_collision_rate_tag.begin(); ite != bcs_collision_rate_tag.end(); ++ite )
    {
        string bc = ite->first;
        int count = ite->second.first;
        double rate = ite->second.second;
        
        if ( bcs_dominant_tags.find(bc) != bcs_dominant_tags.end() )
        {
            string ds = bcs_dominant_tags[bc];
            outf_bt<<bc<<"\t"<<count<<"\t"<<rate<<"\t"<<ds<<endl;
        } else
            outf_bt<<bc<<"\t"<<count<<"\t"<<rate<<"\tCollision"<<endl;
        
    }
    outf_bt.close();

    cout<<"output unknown sample map"<<endl;
    for ( map<string, set<string > >::iterator ite = unknown_sample_map.begin(); ite != unknown_sample_map.end(); ++ite )
    {
        string outfile = ite->first+".unkn_sample_bc.txt";
        ofstream outf( outfile.data() );
        for ( set<string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            outf<<*si<<endl;
        }
        outf.close();
    }
    cout<<"output unknown tag map"<<endl;
    for ( map<string, set<string > >::iterator ite = unknown_tag_map.begin(); ite != unknown_tag_map.end(); ++ite )
    {
        string outfile = ite->first+".unkn_tag_bc.txt";
        ofstream outf( outfile.data() );
        for ( set<string >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            outf<<*si<<endl;
        }
        outf.close();
    }
}

void exit_with_help()
{

	cerr <<"Usage:	Assign_barcode_sample_tag [OPTION1] [VALUE1] [[OPTION2] [VALUE2] ...]" <<endl;
	cerr <<"Options:" <<endl;
	
	cerr <<"-i   input files information in one file: each line recodes one bed file with their meta information; format: [column one] [column two] [column three]" <<endl;
    cerr <<"     column one: one of types:[sa_sample_R2, sa_sample_R1, sa_tag_R2, sa_tag_R1, sc_sample_R2, sc_sample_R1, sc_tag_R2, sc_tag_R1]"<<endl;
    cerr <<"     column two: one of the sample_names or tag_names depending on the file"<<endl;
    cerr <<"     column three: filename"<<endl;
	cerr <<"-r   Collision rate threshold (Optional, default: 0.85)" <<endl;
	cerr <<"-o   output prefix" <<endl;
    cerr <<"-s   for short. Only two types: sa_sample_R2 sa_tag_R2 [sa is the mixture of sa and sc] (default 0 means false; 1 means true)"<<endl;
	
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
    string inputfile = "";
    string outprefix = "out";
    double collision_thr = 0.85;
    int saonly = 0;
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
		
		case 'i':
			inputfile = argv[i];
			break;
        case 'r':
            collision_thr = atof( argv[i] );
            break;
        case 'o':
            outprefix = argv[i];
            break;
        case 's':
            saonly = atoi(argv[i]);
            break;
        
        default:
			exit_with_help();
		}
	}

    if ( inputfile.empty() )
    {
        exit_with_help("Error: inputfile empty!");
    }

    map<string, map<string, map< pair<int, int >, string > > > sample_chr_pos_line_atac;
    map<string, map<string, map<string, map< pair<int, int >, string > > > > sample_tag_chr_pos_line_cuttag;
    map<string, map<string, int > > bc_sample_count;
    map<string, pair<int, double > > bcs_collision_rate_sample;
    map<string, string > bcs_dominant_samples;
    set<string > collision_bcs_sample;
    map<string, map<string, int > > bc_tag_count;
    map<string, pair<int, double > > bcs_collision_rate_tag;
    map<string, string > bcs_dominant_tags;
    set<string> collision_bcs_tag;
    map<string, set<string > > unknown_sample_map;
    map<string, set<string > > unknown_tag_map;

    if ( saonly == 0 )
        readinfiles_main( inputfile, 
            sample_chr_pos_line_atac,
            sample_tag_chr_pos_line_cuttag,
            bc_sample_count,
            bcs_collision_rate_sample,
            bcs_dominant_samples,
            collision_bcs_sample,
            bc_tag_count,
            bcs_collision_rate_tag,
            bcs_dominant_tags,
            collision_bcs_tag,
            unknown_sample_map,
            unknown_tag_map,
            collision_thr );
    else
        readinfiles_main_short( inputfile, 
            sample_chr_pos_line_atac,
            sample_tag_chr_pos_line_cuttag,
            bc_sample_count,
            bcs_collision_rate_sample,
            bcs_dominant_samples,
            collision_bcs_sample,
            bc_tag_count,
            bcs_collision_rate_tag,
            bcs_dominant_tags,
            collision_bcs_tag,
            unknown_sample_map,
            unknown_tag_map,
            collision_thr );
    
    
    output( sample_chr_pos_line_atac, 
        sample_tag_chr_pos_line_cuttag,
        bc_sample_count,
        bcs_collision_rate_sample,
        bcs_dominant_samples,
        collision_bcs_sample,
        bc_tag_count,
        bcs_collision_rate_tag,
        bcs_dominant_tags,
        collision_bcs_tag,
        unknown_sample_map,
        unknown_tag_map,
        outprefix );

    return 1;

}
