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
                } else
                {
                    gene = "unassigned";
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

void readin_cb( string infile, set<string > &bcs, map<string, string > &bc_sample, string prefix )
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
        vector<string > ps2 = parse_string( ps[5], ':');
        string sample = ps2[0];
        bc = prefix+":"+bc;
        bc_sample.insert( make_pair(bc, sample ) );
    }
    inf.close();

}



void readinsam( string infile, set<string > &bcs, map<string, map<string, set<string > > > &cb_gene_umi, 
    map<string, int > &cb_unassigned, string prefix )
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
        if ( i %10000000 == 0 )
            cout<<"processed "<<i<<" lines"<<endl;
        string CB = "";
        string UB = "";
        string gene = "";
        bool pass = parse_line( line, CB, UB,  gene );
        if ( pass )
        {
            if ( bcs.find( CB ) != bcs.end() )
            {
                CB = prefix+":"+CB;
                cb_gene_umi[CB][gene].insert( UB );
            }
        } else
        {
            if ( bcs.find( CB ) != bcs.end() )
            {
                CB = prefix+":"+CB;
                if ( cb_unassigned.find(CB) != cb_unassigned.end() )
                {
                    cb_unassigned[CB] += 1;
                } else
                    cb_unassigned[CB] = 1;

            }
        }
        
    }
    inf.close();
}

void readin_cb_sam_b( string infile, 
    map<string, map<string, set<string > > > &cb_gene_umi, 
    map<string, int > &cb_unassigned, 
    map<string, string > &bc_sample )
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
        string bcfilename = ps[0];
        string samfilename = ps[1];
        string prefix = ps[2];
        set<string > bcs;
        readin_cb( bcfilename, bcs, bc_sample, prefix );
        readinsam( samfilename, bcs, cb_gene_umi, cb_unassigned, prefix );
    }
}

bool check_mediumn( map<string, map<string, int > > &cb_gene_umi, map<string, int > &cb_unassigned, int thr )
{
    bool v = false;
    int total = 0;
    for ( map<string, map<string, int > >::iterator ite = cb_gene_umi.begin(); ite != cb_gene_umi.end(); ++ite )
    {
        for ( map<string, int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            total += si->second;
        }
    }
    for ( map<string, int >::iterator ite = cb_unassigned.begin(); ite != cb_unassigned.end(); ++ite )
    {
        total += ite->second;
    }


    int ave = total / (int)cb_gene_umi.size();
    cout<<"average "<<ave<<" "<<total<<" "<<(int)cb_gene_umi.size()<<" "<<cb_unassigned.size()<<endl;
    
    if ( ave > thr )
    {
        return true;
    } else
        return false;
}


void boost_rna_counts( map<string, map<string, set<string > > > &cb_gene_umi,
    map<string, int > &cb_unassigned,
    set<string > &genes,
    map<string, map<string, int > > &b_cb_gene_umic, 
    map<string, int > &b_cb_unassigned, int thr )
{
    std::srand(std::time(nullptr));
    
    
    for ( map<string, map<string, set<string > > >::iterator ite = cb_gene_umi.begin(); ite != cb_gene_umi.end(); ++ite )
    {
        string cb = ite->first;
        for ( map<string, set<string > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string gene = si->first;
            genes.insert( gene );
            int c = (int)si->second.size();
            b_cb_gene_umic[cb][gene] = c;
        }
    }
    b_cb_unassigned = cb_unassigned;
    int b = 0;
    while ( !check_mediumn( b_cb_gene_umic, b_cb_unassigned, thr ) )
    {
        b+=1;
        if ( b >= 10 )
        {
            cout<<"b=10, break"<<endl;
            break;
        }
        cout<<"boost "<<b<<endl;
        map<string, map<string, int > > tmp_cb_gene_umic;
        map<string, int > tmp_cb_unassigned;
        for ( map<string, map<string, int > >::iterator ite = b_cb_gene_umic.begin(); ite !=  b_cb_gene_umic.end(); ++ite )
        {
            string cb = ite->first;
            for ( set<string >::iterator si = genes.begin(); si != genes.end(); ++si )
            {
                int c = 0;
                if ( ite->second.find( *si ) != ite->second.end() )
                {
                    c = ite->second[*si];
                }
                int add = 0;
                if ( std::rand() % 5 == 0 )
                {
                    add = 1;
                }
                int ca = c*2 + add;
                if ( ca > 0 )
                {
                    tmp_cb_gene_umic[cb][*si] = ca;
                }
            }
            int u = 0;
            if ( b_cb_unassigned.find(cb) != b_cb_unassigned.end() )
            {
                u = b_cb_unassigned[cb] * 2;
            }
            tmp_cb_unassigned[cb] = u;
        }
        b_cb_gene_umic = tmp_cb_gene_umic;
        b_cb_unassigned = tmp_cb_unassigned;

    }


}



void readinepi( string bcfile, string peakfile, string mtxfile, string usfile, string prefix, vector<string > &peak_ve, map<string, map<string, int > > &bc_peak_count, map<string, int > &bc_unassigned )
{
    ifstream inf1( bcfile.data() );
    if( !inf1.good() ){
		std::cout<<"file "<<bcfile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<bcfile<<endl;
	string line;
    vector<string > bc_ve;
    while (!inf1.eof())
	{
        
		getline(inf1,line);
		if ( line.empty() )
			break;
        bc_ve.push_back( line );
    }
    inf1.close();

    ifstream inf2( peakfile.data() );
    if( !inf2.good() ){
		std::cout<<"file "<<peakfile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<peakfile<<endl;
	
    
    while (!inf2.eof())
	{
        
		getline(inf2,line);
		if ( line.empty() )
			break;
        peak_ve.push_back( line );
    }
    inf2.close();

    ifstream inf3( mtxfile.data() );
    if ( !inf3.good() ){
        std::cout<<"file "<<mtxfile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<mtxfile<<endl;
    getline(inf3, line);
    getline(inf3, line );
    while (!inf3.eof())
	{
        
		getline(inf3,line);
		if ( line.empty() )
			break;
        
        vector<string > ps = parse_string( line );
        if ( (int)ps.size() != 3 )
        {
            cout<<"error unexpected line "<<line<<endl; exit(1);
        } 
        int rowi = atoi(ps[0].c_str() );
        int coli = atoi( ps[1].c_str() );
        int c = atoi( ps[2].c_str() );
        if ( rowi > (int)peak_ve.size() )
        {
            cout<<"error unexpected row idx "<<rowi<<" peak_ve "<<peak_ve.size()<<endl; exit(1);
        }
        if ( coli > (int)bc_ve.size() )
        {
            cout<<"error unexpected col idx "<<coli<<" bc_ve "<<bc_ve.size()<<endl; exit(1);
        }
        string peak = peak_ve[rowi-1];
        string bc = bc_ve[coli-1];
        bc = prefix+":"+bc;
        bc_peak_count[bc][peak] = c;
        
    }
    inf3.close();

    ifstream inf4( usfile.data() );
    if ( !inf4.good() ){
        std::cout<<"file "<<usfile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<usfile<<endl;
    while (!inf4.eof())
	{
        
		getline(inf3,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line );
        if ( (int)ps.size() != 2 )
        {
            cout<<"error unexpected line "<<line<<endl; exit(1);
        } 
        string bc = ps[0];
        bc = prefix+":"+bc;
        int c = atoi(ps[1].c_str() );
        bc_unassigned[bc] = c;
    }
    inf4.close();

}



void boost_atac_count( map<string, map<string, int > > &bc_peak_count, 
    map<string, int > &bc_unassigned, 
    vector<string > &peak_ve,
    map<string, map<string, int > > &b_bc_peak_count, 
    map<string, int > &b_bc_unassigned, int thr  )
{
    std::srand(std::time(nullptr));

    b_bc_peak_count = bc_peak_count;
    b_bc_unassigned = bc_unassigned;
    int b = 0;
    while ( !check_mediumn( b_bc_peak_count, b_bc_unassigned, thr ) )
    {
        b+=1;
        if ( b >= 10 )
        {
            cout<<"b=10, break"<<endl;
            break;
        }
        cout<<"boost "<<b<<endl;

        map<string, map<string, int > > tmp_bc_peak_count;
        map<string, int > tmp_bc_unassigned;
        for ( map<string, map<string, int > >::iterator ite = b_bc_peak_count.begin(); ite !=  b_bc_peak_count.end(); ++ite )
        {
            string bc = ite->first;
            for ( size_t i = 0; i < peak_ve.size(); ++i )
            {
                int c = 0;
                string peak = peak_ve[i];
                if ( ite->second.find( peak ) != ite->second.end() )
                {
                    c = ite->second[peak];
                }
                int add = 0;
                if ( std::rand() % 5 == 0 )
                {
                    add = 1;
                }
                int ca = c*2 + add;
                if ( ca > 0 )
                {
                    tmp_bc_peak_count[bc][peak] = ca;
                }
            }
            int u = 0;
            if ( b_bc_unassigned.find(bc) != b_bc_unassigned.end() )
            {
                u = b_bc_unassigned[bc] * 2;
            }
            tmp_bc_unassigned[bc] = u;
        }
        b_bc_peak_count = tmp_bc_peak_count;
        b_bc_unassigned = tmp_bc_unassigned;

    }
}

void readinepi_b( string infile, 
    map<string, vector<string > > &tag_peak_ve, 
    map<string, map<string, map<string, int > > > &tag_b_bc_peak_count,
    map<string, map<string, int > > &tag_b_bc_unassigned )
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
        if ( ps.size() != 7 )
        {
            cout<<"error line "<<line<<endl; exit(1);
        }
        string bcfile = ps[0];
        string peakfile = ps[1];
        string mtxfile = ps[2];
        string usfile = ps[3];
        string prefix = ps[4];
        string tag = ps[5];
        int thr = atoi(ps[6].c_str() );
        vector<string > peak_ve;
        map<string, map<string, int > > bc_peak_count;
        map<string, int > bc_unassigned;
        readinepi( bcfile, peakfile, mtxfile, usfile, prefix, peak_ve, bc_peak_count, bc_unassigned );

        map<string, map<string, int > > b_bc_peak_count;
        map<string, int > b_bc_unassigned;
        boost_atac_count( bc_peak_count, bc_unassigned, peak_ve, b_bc_peak_count, b_bc_unassigned, thr  );
        tag_peak_ve.insert( make_pair(tag, peak_ve ) );
        tag_b_bc_peak_count.insert( make_pair(tag, b_bc_peak_count ) );
        tag_b_bc_unassigned.insert( make_pair(tag, b_bc_unassigned ) );
    }
    inf.close();
}

void boost_cells( map<string, map<string, int > > &b_cb_gene_umic, 
    map<string, int > &b_cb_unassigned, 
    map<string, map<string, map<string, int > > > &tag_b_bc_peak_count,
    map<string, map<string, int > > &tag_b_bc_unassigned,
    map<string, string > &bc_sample,
    int thr )
{
    int original_n = (int)b_cb_gene_umic.size();

    std::srand(std::time(nullptr));
    vector<int > nt;
    for ( int i = 0; i < 100; ++i )
    {
        int w = std::rand() % 20;
        int c = original_n * ( 0.8+ ((w*1.0)/100) );
        
        nt.push_back( c );
    }

   
    std::random_device rd;
    std::mt19937 generator(rd());

    std::vector<std::vector<int>> allResults;
    for ( int i = 0; i < 100; ++i )
    {
        vector<int> numbers;
        for ( int j = 0; j < original_n; j++ )
        {
            numbers.push_back(j);
        }

        std::shuffle(numbers.begin(), numbers.end(), generator);
        vector<int> sn;
        
        for ( int k = 0; k < nt[i]; k++ )
        {
            sn.push_back(numbers[k]);
            
        }
    
        allResults.push_back(sn);
        
    }


    map<string, map<string, int > > add_cb_gene_umic;
    map<string, int > add_cb_unassigned;
    map<string, map<string, map<string, int > > > add_tag_bc_peak_count;
    map<string, map<string, int > > add_tag_bc_unassigned;
    map<string, string > add_bc_sample;
    int lib = 1;
    while ( (int)add_cb_gene_umic.size() + original_n < thr )
    {
        cout<<"iterator "<<lib<<" size"<<(int)add_cb_gene_umic.size() + original_n <<endl;
        stringstream ss;
        ss << lib;
        
        string st = ss.str();
        string prefix = "lib"+st;
        vector<string > bc_ve;
        for ( map<string, map<string, int > >::iterator ite = b_cb_gene_umic.begin(); ite != b_cb_gene_umic.end(); ++ite )
        {
            string bc = ite->first;
            bc_ve.push_back( bc );
        }
        for ( size_t i = 0; i < allResults[lib-1].size(); ++i )
        {
            int idx = allResults[lib-1][i];
            if ( idx >= original_n )
            {
                cout<<"error idx "<<idx<<" larger "<<original_n<<endl;
                exit(1);
            }
            string bc = bc_ve[idx];
            vector<string > ps = parse_string( bc, ':');
            string new_bc = ps[ps.size()-1];
            new_bc = prefix+":"+new_bc;
            map<string, int > n_gene_umic;
            for ( map<string, int >::iterator ite = b_cb_gene_umic[bc].begin(); ite != b_cb_gene_umic[bc].end(); ++ite )
            {
                
                int dif = 5 - (std::rand() % 10);
                int c = ite->second + dif;
                if ( c <= 0 )
                    c = 1;
                n_gene_umic[ite->first] = c;
                
            }
            add_cb_gene_umic[new_bc] = n_gene_umic;

            int dif1 = 5 - (std::rand() % 10);
            int u = b_cb_unassigned[bc] + dif1;
            if ( u <0 )
                u = 0;
            add_cb_unassigned[new_bc] = u;

            for ( map<string, map<string, map<string, int > > >::iterator ite = tag_b_bc_peak_count.begin(); 
                ite != tag_b_bc_peak_count.end(); ++ite )
            {
                map<string, int > n_peak_count;
                for ( map<string, int >::iterator si = ite->second[bc].begin(); si != ite->second[bc].end(); ++si )
                {
                    int dif2 = 5 - (std::rand() % 10);
                    int c = si->second + dif2;
                    if ( c <= 0 )
                        c = 1;
                    n_peak_count[si->first] = c;

                }
                add_tag_bc_peak_count[ite->first][new_bc] = n_peak_count;

                int dif3 = 5 - (std::rand() % 10);
                int u1 = tag_b_bc_unassigned[ite->first][bc] + dif3;
                if ( u1 < 0 )
                    u1 = 0;
                add_tag_bc_unassigned[ite->first][new_bc] = u1;

            }

            string sp = "Unkown";
            if ( bc_sample.find( bc)!= bc_sample.end() )
            {
                sp = bc_sample[bc];
            }
            add_bc_sample[new_bc] = bc;
        }
        lib+=1;
        if ( lib >= 100 )
            break;
    }

    // combine
    cout<<"combine"<<endl;
    for ( map<string, map<string, int > >::iterator ite = add_cb_gene_umic.begin(); ite != add_cb_gene_umic.end(); ++ite )
    {
        b_cb_gene_umic.insert( *ite );
    }
    cout<<" add cb unassigned"<<endl;
  //  b_cb_gene_umic.insert(add_cb_gene_umic.begin(), add_cb_gene_umic.end() );
    for ( map<string, int >::iterator ite = add_cb_unassigned.begin(); ite != add_cb_unassigned.end(); ++ite )
    {
        b_cb_unassigned.insert( *ite );
    } 

 //   b_cb_unassigned.insert( add_cb_unassigned.begin(), add_cb_unassigned.end() );
    cout<<"add bc_peak_count"<<endl;
    for ( map<string, map<string, map<string, int > > >::iterator ite = add_tag_bc_peak_count.begin(); ite != add_tag_bc_peak_count.end(); ++ite )
    {
        string tag = ite->first;
//        tag_b_bc_peak_count[tag].insert(std::make_move_iterator(ite->second.begin()),
//                                std::make_move_iterator(ite->second.end()));
        for ( map<string, map<string, int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            for ( map<string, int >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
                tag_b_bc_peak_count[tag][si->first].insert(*ti);
        } 
        
    }
    cout<<" add tag bc unassigned"<<endl;
    for ( map<string, map<string, int > >::iterator ite = add_tag_bc_unassigned.begin(); ite != add_tag_bc_unassigned.end(); ++ite )
    {   string tag = ite->first;
        for ( map<string, int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
            tag_b_bc_unassigned[tag].insert( *si );
    }
    cout<<" add sample "<<endl;
    for ( map<string, string >::iterator ite = add_bc_sample.begin(); ite != add_bc_sample.end(); ++ite )
    {
        bc_sample.insert( *ite );
    }

}

void output( map<string, map<string, int > > &b_cb_gene_umic, 
    map<string, int > &b_cb_unassigned, 
    set<string > &genes,
    map<string, vector<string > > &tag_peak_ve,
    map<string, map<string, map<string, int > > > &tag_b_bc_peak_count,
    map<string, map<string, int > > &tag_b_bc_unassigned,
    map<string, string > &bc_sample,
    string outprefix )
{
    string outfile_bc = outprefix+".barcodes.txt";
    ofstream outf_bc( outfile_bc.data() );
    vector<string > bc_ve;
    for ( map<string, map<string, int > >::iterator ite = b_cb_gene_umic.begin(); ite != b_cb_gene_umic.end(); ++ite )
    {
        string bc = ite->first;
        bc_ve.push_back( ite->first );
        outf_bc<<bc<<endl;
    }
    outf_bc.close();
    int bc_count = (int)bc_ve.size();

    string outfile_gene = outprefix+".genes.txt";
    ofstream outf_gene( outfile_gene.data() );
    vector<string > gene_ve;
    for ( set<string >::iterator ite = genes.begin(); ite != genes.end(); ++ite )
    {
        string gn = *ite;
        gene_ve.push_back(gn );
        outf_gene<<gn<<endl;
    }
    outf_gene.close();
    int gene_count = (int)gene_ve.size();

    string outfile_rna_mtx = outprefix+".matrix.rna.mtx";
    ofstream outf_rna_mtx( outfile_rna_mtx.data() );
    long umic_total=0;
    for ( size_t i = 0; i < bc_ve.size(); ++i )
    {
        string bc = bc_ve[i];
        for ( size_t j =0; j < gene_ve.size(); ++j )
        {
            string gene = gene_ve[j];
            int c = 0;
            if ( b_cb_gene_umic[bc].find(gene) != b_cb_gene_umic[bc].end() )
                c = b_cb_gene_umic[bc][gene];
            if ( c > 0 )
                umic_total+=c;
        }
    }
    outf_rna_mtx<<"%%MatrixMarket matrix coordinate integer general"<<endl;
    outf_rna_mtx<<gene_count<<"\t"<<bc_count<<"\t"<<umic_total<<endl;
    for ( size_t i = 0; i < bc_ve.size(); ++i )
    {
        string bc = bc_ve[i];
        for ( size_t j =0; j < gene_ve.size(); ++j )
        {
            string gene = gene_ve[j];
            int c = 0;
            if ( b_cb_gene_umic[bc].find(gene) != b_cb_gene_umic[bc].end() )
                c = b_cb_gene_umic[bc][gene];
            if ( c > 0 )
                outf_rna_mtx<<j+1<<"\t"<<i+1<<"\t"<<c<<endl;
        }
    }
    outf_rna_mtx.close();

    string outfile_rna_unassigned = outprefix+".rna.unassigned.txt";
    ofstream outf_rna_una( outfile_rna_unassigned.data() );
    for ( size_t i = 0; i < bc_ve.size(); ++i )
    {
        string bc = bc_ve[i];
        int c = 0;
        if ( b_cb_unassigned.find(bc) != b_cb_unassigned.end() )
        {
            c = b_cb_unassigned[bc];
        }
        outf_rna_una<<i+1<<"\t"<<c<<endl;
    }
    outf_rna_una.close();


 /*   for ( map<string, map<string, map<string, int > > >::iterator ite =  tag_b_bc_peak_count.begin(); ite != tag_b_bc_peak_count.end(); ++ite )
    {
        cout<<ite->first<<endl;
        int i = 0;
        for ( map<string, map<string, int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            cout<<si->first<<endl;
            i +=1;
            if ( i >= 10 )
                break;
            int j = 0;
            for ( map<string, int >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
            {
                cout<<ti->first<<"\t"<<ti->second<<endl;
                j+=1;
                if ( j >= 10 )
                    break;

            }
        }
    } */
    for ( map<string, vector<string > >::iterator ite = tag_peak_ve.begin(); ite != tag_peak_ve.end(); ++ite )
    {
        string tag = ite->first;
        
        int peak_count= (int)ite->second.size();
        string outfile_peak = outprefix+"."+tag+".peaks.txt";
        ofstream outf_pk( outfile_peak.data() );
        for ( size_t i = 0; i < ite->second.size(); ++i )
        {
            outf_pk<<ite->second[i]<<endl;
        }
        outf_pk.close();

        string outfile_peak_mtx = outprefix+"."+tag+".matrix.mtx";
        ofstream outf_pk_mtx( outfile_peak_mtx.data() );
        long count_total = 0;
        for ( size_t i = 0; i < bc_ve.size(); ++i )
        {
            string bc = bc_ve[i];
            for ( size_t j = 0; j < ite->second.size(); ++j )
            {
                string pk = ite->second[j];

                int c= 0;
                if ( tag_b_bc_peak_count[tag][bc].find(pk) != tag_b_bc_peak_count[tag][bc].end() )
                {
                    c = tag_b_bc_peak_count[tag][bc][pk];
                   
                }
                if ( c>0)
                    count_total += c;
            }
        }
        outf_pk_mtx<<"%%MatrixMarket matrix coordinate integer general"<<endl;
        outf_pk_mtx<<peak_count<<"\t"<<bc_count<<"\t"<<count_total<<endl;
        for ( size_t i = 0; i < bc_ve.size(); ++i )
        {
            string bc = bc_ve[i];
            for ( size_t j = 0; j < ite->second.size(); ++j )
            {
                string pk = ite->second[j];
                int c= 0;
                
                if ( tag_b_bc_peak_count[tag][bc].find(pk) != tag_b_bc_peak_count[tag][bc].end() )
                {
                    c = tag_b_bc_peak_count[tag][bc][pk];
                }
                if ( c>0)
                    outf_pk_mtx<<j+1<<"\t"<<i+1<<"\t"<<c<<endl;
            }
        }
        outf_pk_mtx.close();

        string outfile_peak_unassigned = outprefix+"."+tag+".unassigned.txt";
        ofstream outf_pk_una( outfile_peak_unassigned.data() );
        for ( size_t i = 0; i < bc_ve.size(); ++i )
        {
            string bc = bc_ve[i];
            int c= 0;
            if ( tag_b_bc_unassigned[tag].find(bc) != tag_b_bc_unassigned[tag].end() )
            {
                c = tag_b_bc_unassigned[tag][bc];
            }
            if ( c > 0)
                outf_pk_una<<i+1<<"\t"<<c<<endl;
        }
        outf_pk_una.close();

    }

    string outfile_sample = outprefix+".sample.txt";
    ofstream outf_sp( outfile_sample.data() );
    for ( size_t i = 0; i < bc_ve.size(); ++i )
    {
        string bc = bc_ve[i];
        string sp = "unknown";
        if ( bc_sample.find( bc) == bc_sample.end() )
        {
            sp = bc_sample[bc];
        }
        outf_sp<<i+1<<"\t"<<sp<<endl;
    }
    outf_sp.close();
}

int main( int argc, char* argv[] )
{
    if ( argc  == 1 )
    {
        cout<<"Integrate all the materials"<<endl;
        cout<<"Usage: prog file_bc_sam[b] file_epi[b] outprefix"<<endl;
        exit(1);
    }

    int thr_count = 5000;
    int thr_cell = 50000;

    string infile_bc_sam = argv[1];
    string infile_epi = argv[2];
    string outprefix = argv[3];

    map<string, map<string, set<string > > > cb_gene_umi;
    map<string, int > cb_unassigned;
    map<string, string > bc_sample;
    readin_cb_sam_b( infile_bc_sam, cb_gene_umi, cb_unassigned, bc_sample );

    cout<<"boost"<<endl;
    set<string > genes;
    map<string, map<string, int > > b_cb_gene_umic;
    map<string, int > b_cb_unassigned;
    boost_rna_counts( cb_gene_umi, cb_unassigned, genes, b_cb_gene_umic,  b_cb_unassigned, thr_count );

    cout<<"readin epi"<<endl;
    map<string, vector<string > > tag_peak_ve;
    map<string, map<string, map<string, int > > > tag_b_bc_peak_count;
    map<string, map<string, int > > tag_b_bc_unassigned;
    readinepi_b( infile_epi, 
        tag_peak_ve, 
        tag_b_bc_peak_count,
        tag_b_bc_unassigned );
  
    cout<<"boost cells"<<endl;
    boost_cells( b_cb_gene_umic, b_cb_unassigned, tag_b_bc_peak_count, tag_b_bc_unassigned, bc_sample, thr_cell );

    cout<<"Output"<<endl;
    output( b_cb_gene_umic, 
        b_cb_unassigned, 
        genes,
        tag_peak_ve,
        tag_b_bc_peak_count,
        tag_b_bc_unassigned,
        bc_sample,
        outprefix );

    return 1;
}




