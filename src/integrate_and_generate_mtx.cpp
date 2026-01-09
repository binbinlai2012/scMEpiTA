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
#include <omp.h>

using namespace std;

string inttostr( int i )
{
    stringstream ss;
    ss << i;
    
    string st = ss.str();
    return st;
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

vector<string > parse_string( const string & instr, char spl )
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

bool parse_filter_bcline( string line, size_t fsize, map<size_t, string >& idx_tag,
    string &bc, string &sample, string &tag, double colrate_thr, int tag_nthr )
{
    bool pas = false;
    vector<string > ps = parse_string( line );
    if ( ps.size() != fsize )
    {
        cout<<"unexpected line "<<line<<endl; exit(1);
    }
    double colrate = atof( ps[4].c_str() );
    
    if ( colrate < colrate_thr )
        return false;
    vector< string > ps2 = parse_string( ps[5], ':');
    bc = ps[0];
    sample = ps2[0];
    size_t max_tag_idx = 0;
    int max_tag_n = 0;
    int total_tag_n = 0;
    for ( size_t i = 6; i < fsize; ++i )
    {
        int tag_n = atoi( ps[i].c_str() );
        if ( tag_n > max_tag_n )
        {
            max_tag_n = tag_n;
            max_tag_idx = i;
        }
        total_tag_n += tag_n;
    }
    if ( max_tag_n > tag_nthr )
    {
        double r = max_tag_n*1.0/total_tag_n;
        
        if ( r < colrate_thr )
            return false;
        else
        {
            tag = idx_tag[max_tag_idx];
        }
    }
    return true;

}

void readin_and_filter_cb( string infile, 
    set<string > &bcs, 
    map<string, string > &bc_sample, 
    map<string, string > &bc_tag, 
    double colrate_thr, 
    int tag_nthr )
{
    ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
    getline(inf, line);

    map<size_t, string > idx_tag;
    vector<string > ps = parse_string( line );
    size_t fsize = ps.size();
    for ( size_t i = 6; i < ps.size(); ++i )
        idx_tag.insert( make_pair(i, ps[i]) );

    string outfile_ab = infile+".failed.txt";
    ofstream outf( outfile_ab.data() );
    outf<<line<<endl;

    while(!inf.eof())
	{
        
		getline(inf,line);
		if ( line.empty() )
			break;
        string bc = "";
        string sample = "";
        string tag = "";
        bool pas = parse_filter_bcline( line, fsize, idx_tag, 
            bc, sample, tag, colrate_thr, tag_nthr );
        if ( pas )
        {
            bcs.insert( bc );
            bc_sample[bc] = sample;
            if ( tag != "" )
                bc_tag[bc] = tag;
        } else
            outf<<line<<endl;
    }
    inf.close();
    outf.close();

}

void readinsam( string infile, set<string > &bcs, map<string, map<string, set<string > > > &cb_gene_umi, 
    map<string, int > &cb_unassigned )
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
                
                cb_gene_umi[CB][gene].insert( UB );
            }
        } else
        {
            if ( bcs.find( CB ) != bcs.end() )
            {
               
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


void readinsam(string infile, set<string> &bcs, map<string, map<string, set<string>>> &cb_gene_umi, 
    map<string, int> &cb_unassigned, int threads = 0)
{
    ifstream inf(infile.data());
    if(!inf.good()){
        std::cout<<"file "<<infile<<" not found!"<<std::endl;
        exit(1);
    }
    cout<<"read file "<<infile<<endl;
    
    // 设置线程数
    if(threads > 0) {
        omp_set_num_threads(threads);
        cout << "Using " << threads << " threads" << endl;
    } else {
        int max_available = omp_get_max_threads();
        int suggested = max_available;
        if(max_available > 8) {
            suggested = 8;
        }
        omp_set_num_threads(suggested);
        cout << "Using " << suggested << " threads (auto-detected from " << max_available << " available)" << endl;
    }
    
    // 读取所有行到内存
    vector<string> lines;
    string line;
    while(getline(inf, line)) {
        if(!line.empty() && line[0] != '@') {
            lines.push_back(line);
        }
    }
    inf.close();
    
    cout<<"Total lines to process: "<<lines.size()<<endl;
    
    // 获取线程数
    int num_threads = omp_get_max_threads();
    
    // 线程本地存储
    vector<map<string, map<string, set<string>>>> thread_cb_gene_umi(num_threads);
    vector<map<string, int>> thread_cb_unassigned(num_threads);
    
    // 记录开始时间
    double start_time = omp_get_wtime();
    
    // 进度计数器
    size_t processed_count = 0;
    size_t report_interval = max((size_t)1, lines.size() / 20); // 报告20次
    
    // 并行处理
    #pragma omp parallel reduction(+:processed_count)
    {
        int thread_id = omp_get_thread_num();
        
        #pragma omp for schedule(dynamic, 10000)
        for(size_t i = 0; i < lines.size(); ++i) {
            // 处理逻辑
            const string& line = lines[i];
            string CB = "";
            string UB = "";
            string gene = "";
            
            bool pass = parse_line(line, CB, UB, gene);
            if(pass) {
                if(bcs.find(CB) != bcs.end()) {
                    thread_cb_gene_umi[thread_id][CB][gene].insert(UB);
                }
            } else {
                if(bcs.find(CB) != bcs.end()) {
                    thread_cb_unassigned[thread_id][CB] += 1;
                }
            }
            
            // 更新计数器
            processed_count++;
            
            // 进度报告
            if(processed_count % report_interval == 0) {
                #pragma omp critical
                {
                    float percent = 100.0 * processed_count / lines.size();
                    cout << fixed << setprecision(1) << "Progress: " << percent 
                         << "% (" << processed_count << "/" << lines.size() << ")" << endl;
                }
            }
        }
    }
    
    // 记录并行处理结束时间
    double parallel_end_time = omp_get_wtime();
    cout << "Parallel processing time: " << (parallel_end_time - start_time) << " seconds" << endl;
    
    // 合并结果
    double merge_start_time = omp_get_wtime();
    
    for(int t = 0; t < num_threads; ++t) {
        for(const auto& cb_entry : thread_cb_gene_umi[t]) {
            const string& cb = cb_entry.first;
            for(const auto& gene_entry : cb_entry.second) {
                const string& gene = gene_entry.first;
                cb_gene_umi[cb][gene].insert(gene_entry.second.begin(), gene_entry.second.end());
            }
        }
        
        for(const auto& unassigned_entry : thread_cb_unassigned[t]) {
            cb_unassigned[unassigned_entry.first] += unassigned_entry.second;
        }
    }
    
    double merge_end_time = omp_get_wtime();
    cout << "Result merging time: " << (merge_end_time - merge_start_time) << " seconds" << endl;
    cout << "Total time: " << (merge_end_time - start_time) << " seconds" << endl;
    cout << "Processed " << processed_count << " lines in total" << endl;
}

void readin_cb_sam_b( string infile, 
    set<string > &bcs, 
    map<string, map<string, set<string > > > &cb_gene_umi, 
    map<string, int > &cb_unassigned, 
    map<string, string > &bc_sample,
    map<string, string > &bc_tag,
    double colrate_thr, 
    int tag_nthr,
    int threads )
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
        readin_and_filter_cb( bcfilename, bcs, bc_sample, bc_tag, colrate_thr, tag_nthr );

        readinsam( samfilename, bcs, cb_gene_umi, cb_unassigned, threads);
    }

	inf.close();
}

void readinpeaks( string infile, map<string, set<pair<int, int > > > &peaks, map<string, int> &peak_id, vector<string > &peak_ve )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
    int idx = 0;
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
        string p = chr+"-"+inttostr(st)+"-"+inttostr(ed);
        peak_ve.push_back(p);
        idx+=1;
        peak_id[p] = idx;
	}
	inf.close();
}
/*
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
} */

string cal_overlap(const string &chr, int st, int ed, const map<string, set<pair<int, int>>> &peaks)
{
    string p = "None";
    bool fd = false;
    
    // 使用find方法查找chr，避免使用operator[]
    auto it = peaks.find(chr);
    if (it == peaks.end()) {
        return p;  // 如果没有找到染色体，直接返回"None"
    }
    
    // 如果找到了，使用it->second来获取该染色体对应的set
    const set<pair<int, int>> &peak_set = it->second;
    for (set<pair<int, int>>::const_iterator ite = peak_set.begin(); ite != peak_set.end(); ++ite)
    {
        if (ite->second < st)
            continue;
        if (ite->first > ed)
            break;
        fd = true;
        p = chr + "-" + inttostr(ite->first) + "-" + inttostr(ite->second);
        break;
    }
    return p;
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
void assign_read_to_peaks( string infile, 
    string tag,
	map<string, string > &bc_tag, 
	map<string, set<pair<int, int > > > &peaks, 
    map<string, int > &peak_idx,
	map<string, map<int, int > > &bc_peakid_count,
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
        if ( tag != "atac" )
        {
            if ( bc_tag.find(bc) == bc_tag.end() )
                continue;
            else
            {
                
                if ( bc_tag[bc] != tag  )
                    continue;
            }
        }
        

        string chr = ps[0];
        
        int st = atoi(ps[1].c_str() );
		int ed = atoi(ps[2].c_str() );

        string ovlp = cal_overlap( chr, st, ed, peaks );
        if ( ovlp != "None")
        {
            if ( peak_idx.find( ovlp ) != peak_idx.end() )
            {
                int idx = peak_idx[ovlp];
                if ( bc_peakid_count[bc].find( idx ) != bc_peakid_count[bc].end() )
                {
                    bc_peakid_count[bc][idx] += 1;
                } else
                {
                    bc_peakid_count[bc][idx] = 1;
                }
            } else
            {
                cout<<"error cannot find peak in peak_idx "<<ovlp<<endl; exit(1);
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

} */

void assign_read_to_peaks( string infile, 
    string tag,
    set<string > &bcs, 
    map<string, string > &bc_tag, 
    map<string, set<pair<int, int > > > &peaks, 
    map<string, int > &peak_idx,
    map<string, map<int, int > > &bc_peakid_count,
    map<string, int > &unmapped,
    int threads=0 )
{
    ifstream inf(infile.data());
    if(!inf.good()){
        std::cout<<"file "<<infile<<" not found!"<<std::endl;
        exit(1);
    }
    cout<<"read file "<<infile<<endl;
    
    // 设置线程数
    if(threads > 0) {
        omp_set_num_threads(threads);
        cout << "Using " << threads << " threads" << endl;
    } else {
        // 自动检测，但限制最大线程数避免过度消耗资源
        int max_available = omp_get_max_threads();
        // 通常建议使用CPU核心数，可以自行调整
        int suggested = max_available;
        if(max_available > 8) {
            suggested = 8; // 限制最大线程数，避免过多线程导致上下文切换开销
        }
        omp_set_num_threads(suggested);
        cout << "Using " << suggested << " threads (auto-detected from " << max_available << " available)" << endl;
    }
    
    // 读取所有行到内存
    vector<string> lines;
    string line;
    while(getline(inf, line)) {
        if(!line.empty()) {
            lines.push_back(line);
        }
    }
    inf.close();
    
    cout<<"Total lines: "<<lines.size()<<endl;
    
    // 获取当前设置的线程数
    int num_threads = omp_get_max_threads();
    
    // 线程本地存储
    vector<map<string, map<int, int>>> thread_results(num_threads);
    vector<map<string, int>> thread_unmapped(num_threads);
    
    // 创建只读引用
    const auto& bc_tag_ref = bc_tag;
    const auto& peaks_ref = peaks;
    const auto& peak_idx_ref = peak_idx;
    const auto& bcs_ref = bcs;
    
    // 记录开始时间
    double start_time = omp_get_wtime();
    
    // 并行处理
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        
        #pragma omp for schedule(dynamic, 1000)
        for(size_t i = 0; i < lines.size(); ++i) {
            vector<string> ps = parse_string(lines[i], '\t');
            
            if(ps.size() < 4) continue;
            
            string bc = parse_bc(ps[3]);

            auto bc_it0 = bcs_ref.find(bc);
            if ( bc_it0 == bcs_ref.end() )
                continue;
            
            if ( tag != "atac" )
            {
                auto bc_it = bc_tag_ref.find(bc);
                if(bc_it == bc_tag_ref.end() || bc_it->second != tag) {
                    continue;
                }
            }
            
            
            string chr = ps[0];
            int st_val = atoi(ps[1].c_str());
            int ed_val = atoi(ps[2].c_str());
            
            string ovlp = cal_overlap(chr, st_val, ed_val, peaks_ref);
            if(ovlp != "None") {
                auto idx_it = peak_idx_ref.find(ovlp);
                if(idx_it != peak_idx_ref.end()) {
                    thread_results[thread_id][bc][idx_it->second] += 1;
                } else {
                    #pragma omp critical
                    {
                        cerr << "Warning: cannot find peak in peak_idx " << ovlp << endl;
                    }
                }
            } else {
                thread_unmapped[thread_id][bc] += 1;
            }
        }
    }
    
    // 记录结束时间
    double end_time = omp_get_wtime();
    cout << "Parallel processing time: " << (end_time - start_time) << " seconds" << endl;
    
    // 合并结果
    start_time = omp_get_wtime();
    for(int t = 0; t < num_threads; ++t) {
        for(const auto& bc_entry : thread_results[t]) {
            const string& bc = bc_entry.first;
            for(const auto& peak_entry : bc_entry.second) {
                bc_peakid_count[bc][peak_entry.first] += peak_entry.second;
            }
        }
        
        for(const auto& unmapped_entry : thread_unmapped[t]) {
            unmapped[unmapped_entry.first] += unmapped_entry.second;
        }
    }
    end_time = omp_get_wtime();
    cout << "Result merging time: " << (end_time - start_time) << " seconds" << endl;
}

void readinepi( set<string > &bcs,
    map<string, string > &bc_tag,
    string peakfile,    // tag filename
    string bedfile,     // tag filename
    map<string, vector<string > > &tag_peak_ve, 
    map<string, map<string, map<int, int > > > &tag_bc_peak_count, 
    map<string, map<string, int > > &tag_bc_unassigned,
    int threads )
{
    // read in peak files
    ifstream inf1( peakfile.data() );
    if( !inf1.good() ){
		std::cout<<"file "<<peakfile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<peakfile<<endl;
	string line;
    map<string, map<string, set<pair<int, int > > > > tag_peaks;
    map<string, map<string, int > > tag_peak_idx;
    while (!inf1.eof())
	{
        
		getline(inf1,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line );
        string tag = ps[0];
        string filename = ps[1];
        map<string, set<pair<int, int > > > peaks;
        map<string, int > peak_idx;
        vector<string > peakve;
        readinpeaks( filename, peaks, peak_idx, peakve );
        tag_peaks[tag] = peaks;
        tag_peak_idx[tag] = peak_idx;
        tag_peak_ve[tag] = peakve;
    }
    inf1.close();

    ifstream inf2( bedfile.data() );
    if( !inf2.good() ){
		std::cout<<"file "<<bedfile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<bedfile<<endl;
	
    map<string, string > tag_bedfilename;
    while (!inf2.eof())
	{
        
		getline(inf2,line);
		if ( line.empty() )
			break;
        vector<string > ps = parse_string( line );
        string tag = ps[0];
        string filename = ps[1];
        if ( tag_peaks.find(tag) == tag_peaks.end() )
        {
            cout<<"error cannot find tag in tag_peaks "<<tag<<endl; exit(1);
            
        }
        if ( tag_peak_idx.find( tag ) == tag_peak_idx.end() )
        {
            cout<<"error cannot find tag in tag_peak_idx "<<tag<<endl; exit(1);
        }

        map<string, map<int, int > > bc_peakid_count;
        map<string, int > unmapped;
        assign_read_to_peaks( filename, tag, bcs, bc_tag,  tag_peaks[tag], tag_peak_idx[tag], bc_peakid_count, unmapped, threads );

        tag_bc_peak_count[tag] = bc_peakid_count;
        tag_bc_unassigned[tag] = unmapped;
    }
    inf2.close();

}

void output( string outprefix,
    set<string > &bcs, 
    map<string, string > &bc_sample,
    map<string, map<string, set<string > > > &cb_gene_umi, 
    map<string, int > &cb_unassigned, 
    map<string, vector<string > > &tag_peak_ve, 
    map<string, map<string, map<int, int > > > &tag_bc_peak_count, 
    map<string, map<string, int > > &tag_bc_unassigned )
{
    if ( outprefix[outprefix.size()-1] != '/' )
        outprefix+='.';
    map<string, int > bc_idx;
    int idx = 0;
    string outfile_bc = outprefix+"barcodes.txt";
    string outfile_bcsample = outprefix+"bc_sample.txt";
    ofstream outf_bc(outfile_bc.data() );
    ofstream outf_bcsamp( outfile_bcsample.data() );
    for ( set<string >::iterator ite = bcs.begin(); ite != bcs.end(); ++ite )
    {
        idx+= 1;
        bc_idx.insert( make_pair(*ite, idx) );
        outf_bc << *ite <<endl;
        outf_bcsamp<<idx<<"\t"<<bc_sample[*ite]<<endl;
    }
    outf_bc.close();   
    outf_bcsamp.close();
   

    set<string > genes;
    for ( map<string, map<string, set<string > > >::iterator ite = cb_gene_umi.begin(); ite != cb_gene_umi.end(); ++ite)
    {
     
        for ( map<string, set<string > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            genes.insert( si->first );
        }
    }

    map<string, int > gene_idx;
    idx = 0;
    string outfile_gene = outprefix+"genes.txt";
    ofstream outf_gene( outfile_gene.data() );
    for ( set<string >::iterator ite = genes.begin(); ite != genes.end(); ++ite )
    {
        idx+=1;
        gene_idx.insert(make_pair(*ite, idx));
        outf_gene<<*ite<<endl;
    }
    outf_gene.close();

    string outfile_summary = outprefix+"summary.txt";
    ofstream outf_summary( outfile_summary.data() );
    outf_summary<<"CB\tType\tCount"<<endl; 

    cout<<"gene mtx"<<endl;
    string outfile_genemtx = outprefix+"gene.mtx";
//    string outfile_genemtx_h = outfile_genemtx+".header";
    ofstream outf_genemtx( outfile_genemtx.data() );
    outf_genemtx<<"%%MatrixMarket matrix coordinate integer general"<<endl;
//    ofstream outf_genemtx_h( outfile_genemtx_h.data() );
//    outf_genemtx_h<<"%%MatrixMarket matrix coordinate integer general"<<endl;
    long total_c_umi = 0;
    vector<string > genemtx_ve;
    for ( map<string, map<string, set<string > > >::iterator ite = cb_gene_umi.begin(); ite != cb_gene_umi.end(); ++ite)
    {
        string bc = ite->first;
        int idxb = bc_idx[bc];
        int bc_umic = 0;
        for ( map<string, set<string > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
        {
            string gene = si->first;
            int idxg = gene_idx[gene];
            int umic = (int)si->second.size();
        //    outf_genemtx<<idxg<<"\t"<<idxb<<"\t"<<umic<<endl;
            genemtx_ve.push_back(inttostr(idxg)+"\t"+inttostr(idxb)+"\t"+inttostr(umic) );
            bc_umic+=umic;
        }    
        outf_summary<<idxb<<"\tUMIC_ingene\t"<<bc_umic<<endl;
        outf_summary<<idxb<<"\tGeneC\t"<<ite->second.size()<<endl;
        total_c_umi += (int)ite->second.size();

        int uns = 0;
        if (cb_unassigned.find( bc ) != cb_unassigned.end() )
            uns += cb_unassigned[bc];
        outf_summary<<idxb<<"\tUMIC_total\t"<<uns+bc_umic<<endl;
    }
    outf_genemtx<<genes.size()<<"\t"<<bcs.size()<<"\t"<<total_c_umi<<endl;
    for ( size_t i = 0; i < genemtx_ve.size(); ++i )
        outf_genemtx<<genemtx_ve[i]<<endl;
    vector<string>().swap(genemtx_ve);
    outf_genemtx.close();
 //   outf_genemtx_h<<genes.size()<<"\t"<<bcs.size()<<"\t"<<total_c<<endl;

    set<string > tags;
    for ( map<string, vector<string > >::iterator ite = tag_peak_ve.begin(); ite != tag_peak_ve.end(); ++ite )
    {
        tags.insert( ite->first);
    }

    cout<<"tag"<<endl;
    
    for ( set<string >::iterator ii = tags.begin(); ii != tags.end(); ++ii )
    {
        string tag = *ii;
        cout<<tag<<endl;
        string outfile_tag_peaks = outprefix+tag+".peaks.txt";
        ofstream outf_tag_peaks( outfile_tag_peaks.data() );
        for ( vector<string >::iterator ite = tag_peak_ve[tag].begin(); ite != tag_peak_ve[tag].end(); ++ite )
        {
            outf_tag_peaks<<*ite<<endl;
        }
        outf_tag_peaks.close();

        string outfile_peak_count = outprefix+tag+".peaks.mtx";
//        string outfile_peak_count_h = outfile_peak_count+".header";
        ofstream outf_peak_count( outfile_peak_count.data() );
//        ofstream outf_peak_count_h( outfile_peak_count_h.data() );
        long total_c = 0;
        set<int > idx_bcs;
        map<int, int > bc_countinpeak;
        map<int, int > bc_countoutpeak;
        vector<string > tag_mtx_ve;
        for ( map<string, map<int, int > >::iterator ite = tag_bc_peak_count[tag].begin(); ite != tag_bc_peak_count[tag].end(); ++ite )
        {
            string bc = ite->first;
            if ( bc_idx.find(bc) == bc_idx.end() )
            {
                cout<<"error cannot find bc "<<bc<<endl; exit(1);
            }
            int idxb = bc_idx[bc];
            idx_bcs.insert( idxb );
            
            int tc = 0;
            for ( map<int, int >::iterator si = ite->second.begin(); si != ite->second.end();++si )
            {
                int idxp = si->first;
                int c = si->second;
        //        outf_peak_count<<idxp<<"\t"<<idxb<<"\t"<<c<<endl;
                tag_mtx_ve.push_back( inttostr(idxp)+"\t"+inttostr(idxb)+"\t"+inttostr(c) );
                total_c += 1;
                tc += c;
            }
            bc_countinpeak[idxb] = tc;
        }
        outf_peak_count<<"%%MatrixMarket matrix coordinate integer general"<<endl;
        outf_peak_count<<tag_peak_ve[tag].size()<<"\t"<<bcs.size()<<"\t"<<total_c<<endl;
        for ( size_t i = 0; i < tag_mtx_ve.size(); ++i )
            outf_peak_count<<tag_mtx_ve[i]<<endl;
        outf_peak_count.close();
     //   outf_peak_count_h<<"%%MatrixMarket matrix coordinate integer general"<<endl;
     //   outf_peak_count_h<<tag_peak_ve[tag].size()<<"\t"<<bcs.size()<<"\t"<<total_c<<endl;
     //   outf_peak_count_h.close();

        string outfile_unas = outprefix+tag+".out_of_peaks.txt";
        ofstream outf_unas( outfile_unas.data() );
        for ( map<string, int >::iterator ite = tag_bc_unassigned[tag].begin(); ite != tag_bc_unassigned[tag].end(); ++ite )
        {
            string bc = ite->first;
            if ( bc_idx.find(bc) == bc_idx.end() )
            {
                cout<<"error cannot find bc "<<bc<<endl; exit(1);
            }
            int idxb = bc_idx[bc];
            outf_unas<<idxb<<"\t"<<ite->second<<endl;
            bc_countoutpeak[idxb] = ite->second;
        }
        outf_unas.close();

        for ( set<int >::iterator ite = idx_bcs.begin(); ite != idx_bcs.end(); ++ite )
        {
            
            int cinp = 0;
            if ( bc_countinpeak.find( *ite ) != bc_countinpeak.end() )
                cinp = bc_countinpeak[*ite];
            int coutp = 0;
            if ( bc_countoutpeak.find( *ite ) != bc_countoutpeak.end() )
                coutp = bc_countoutpeak[*ite];
            outf_summary<<*ite<<"\t"<<tag<<"_inpeakC\t"<<cinp<<endl;
            outf_summary<<*ite<<"\t"<<tag<<"_totalC\t"<<coutp+cinp<<endl;
            
        }
    }
    outf_summary.close();


}

int main( int argc, char* argv[] )
{
    if ( argc == 1 )
    {
        cout<<"Integrate multi-omics and generate matrix"<<endl;
        cout<<"Usage: Prog bc_sam peakfile bedfile outprefix threads"<<endl;
        exit(1);
    }

    string inbcsamfile = argv[1];
    string peakfile = argv[2];
    string bedfile= argv[3];
    string prefix = argv[4];
    int threads = 0;
    if ( argc == 6 )
        threads = atoi(argv[5]);

    double colrate_thr = 0.85;
    int tag_nthr = 200;

    cout<<" read in cb sam files"<<endl;
    set<string > bcs;
    map<string, map<string, set<string > > > cb_gene_umi;
    map<string, int > cb_unassigned;
    map<string, string > bc_sample;
    map<string, string > bc_tag;
    readin_cb_sam_b( inbcsamfile, bcs, cb_gene_umi, cb_unassigned, bc_sample, bc_tag, colrate_thr, tag_nthr, threads );

    cout<<" read in peak bedfile"<<endl;
    map<string, vector<string > > tag_peak_ve;
    map<string, map<string, map<int, int > > > tag_bc_peak_count;
    map<string, map<string, int > > tag_bc_unassigned;
    readinepi( bcs, bc_tag, peakfile,  bedfile, tag_peak_ve, tag_bc_peak_count,  tag_bc_unassigned, threads );

    cout<<"output"<<endl;
    output( prefix, bcs, bc_sample, cb_gene_umi, cb_unassigned, tag_peak_ve, tag_bc_peak_count, tag_bc_unassigned );

    return 1;

}









