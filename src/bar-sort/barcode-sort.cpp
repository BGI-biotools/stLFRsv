#include <stdio.h>
#include <stdlib.h>
#include <htslib/sam.h>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <pthread.h>
#include <signal.h>
#include <unistd.h>
#include <iomanip>
#include <semaphore.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <assert.h>
#include <vector>
#include <time.h>
#include <regex>
#include <cmath>

#include "barcode-sort.h"

// using boost::lexical_cast;
// using boost::bad_lexical_cast;
using namespace std;


int main(int argc,char* argv[]){
	if(argc!=10){
		cerr <<"Usage:"<<argv[0]<<" [LFR_bam_file] [thread_num] [out_directory] [out_prefix] [bar-threshold] [seg-threshold] [gap-size] [read-pair size] [filter-chr(Y/N)]"<<endl;
		cerr <<"Sort the barcoded LFR bam"<<endl;
		cerr <<"Version 0.3"<<endl;
		exit(1);
	}
	string bam_file=argv[1];
	int thread=atoi(argv[2]);
	string dir=argv[3];
	string prefix=argv[4];
	uint32_t bar_th=atoi(argv[5]);
	uint32_t seg_th=atoi(argv[6]);
	uint32_t gap_th=atoi(argv[7]);
	uint32_t is_th=atoi(argv[8]);
	char flag=*argv[9];
	
	//read chromesome
	samFile *in = sam_open(bam_file.c_str(), "r");
	bam_hdr_t *header = sam_hdr_read(in); 
	vector<chr> chrom;
	int tid_count=header->n_targets;
	if(tid_count <1){
		cerr<<"bam header error,chromesome name is missing!!"<<endl;
		exit(1);
	}
	for(int i=0;i<tid_count;i++){
		chr tmp_chrom;
		tmp_chrom.tid_len=(header->target_len)[i];
		tmp_chrom.index_name=i;
		chrom.push_back(tmp_chrom);
	}
	
	sort(chrom.begin(),chrom.end(),comp_chrom);
	uint32_t max_len=chrom.at(0).tid_len * 1.1;
	
	//add the num index associated with chromesome name
	// for(int i=0;i < chrom.size();i++){
		// chrom.at(i).index_name=(uint8_t) i;	
	// }
	
	vector<vector <chr> > sets (chrom.size());
	sets.at(0).push_back(chrom.at(0));
	uint32_t m=1;
	uint32_t n=1;
	while(n<chrom.size()){
		//add the front one;
		m++;
		uint32_t add_len=0;
		sets.at(m-1).push_back(chrom.at(n));
		add_len+=chrom.at(n).tid_len;
		
		//if possible,add the last one
		while(n<(chrom.size()-1)){
			add_len+=chrom.back().tid_len;
			if(add_len <= max_len){
				sets.at(m-1).push_back(chrom.back());
				chrom.pop_back();
			}else{
				break;
			}
		}
		n++;
	}
	// bam_hdr_destroy(header);
	sam_close(in); 
	//read chromesome
	
	cout<<"Get the chromesome info...Done"<<endl;
	//prepare multi threads
	vector<pthread_t> tid(thread);
	vector<int> tid_id(thread);
	for(int i=0;i<thread;i++){
		tid_id.at(i)=i;
	}
	sem_init(&block_sem, 0, thread);
	pthread_mutex_init(&mutex1,NULL);
	pthread_attr_t attr; 
    pthread_attr_init (&attr);
	pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_DETACHED);
	vector<sp_arg> arg_sp(thread);
	
	//split the file and get the info
	int m_start = time(NULL);
	int total_time=m_start;
	cout<<"Split and get the proper pair from the file..."<<endl;
	cout<<"Running in "<< sets.size()<<" groups with "<<thread<<" threads..."<<endl;
	vector<mem_block> mem_cache(sets.size());
	for(uint32_t i=0;i<sets.size();i++){
		sem_wait(&block_sem);
		usleep(5000);
		int j;
		//init the mem cache
		mem_cache.at(i).count=0;
		mem_cache.at(i).size=INIT_SPLIT;
		mem_cache.at(i).l_mem=(read_pair *) malloc(INIT_SPLIT * sizeof(read_pair));
		//
		pthread_mutex_lock(&mutex1);
		j=tid_id.back();
		tid_id.pop_back();
		pthread_mutex_unlock(&mutex1);
		arg_sp.at(j).th_index=j;
		arg_sp.at(j).bam_file=bam_file;
		arg_sp.at(j).chrom=sets.data()+i;
		arg_sp.at(j).mem=mem_cache.data()+i;
		arg_sp.at(j).tid_id=&tid_id;
		arg_sp.at(j).flag=flag;
		arg_sp.at(j).is_th=is_th;
		int res=pthread_create(&tid.at(j),&attr,split_process,(void*)(arg_sp.data()+j));
		usleep(5000);
		if(res!=0){
			cerr<<"pthread_create error:error code="<<res<<endl;
			exit(1);
		}
	}
	
	while(1){
		int value;
		sem_getvalue(&block_sem,&value);
		if(value == thread){
			break;
		}
	}
	cout<<"Split and get the proper pair from the file...Done"<<endl;
	int m_end=time(NULL);
	cout<< "Elapsed time " << m_end - m_start << " seconds" << endl;
	m_start=m_end;
	
	cout<<"Merge and Sort the pair alignment info..."<<endl;
	
	uint64_t total_count=0;
	for(uint32_t i=0;i<sets.size();i++){
		total_count+=mem_cache.at(i).count;
	}
	
	mem_block total_pair;
	total_pair.count=total_count;
	total_pair.size=total_count;
	total_pair.l_mem=(read_pair *) malloc(total_count * sizeof(read_pair));
	
	//copy the mem_cache
	uint64_t mark=0;
	for(uint32_t i=0;i<sets.size();i++){
		memcpy(total_pair.l_mem+mark,mem_cache.at(i).l_mem,mem_cache.at(i).count * sizeof(read_pair));
		mark+=mem_cache.at(i).count;
		free(mem_cache.at(i).l_mem);
	}
	vector<mem_block>().swap(mem_cache);
	//sort
	sort(total_pair.l_mem,total_pair.l_mem+total_count,comp_read);
	cout<<"Merge and Sort the pair alignment info...Done"<<endl;
	m_end=time(NULL);
	cout<< "Elapsed time " << m_end - m_start << " seconds" << endl;
	m_start=m_end;
	
	cout<<"Get the LFR read segments..."<<endl;
	//get the seg
	uint64_t split_count=(total_count/thread)+1;
	mark=0;
	vector<seg_arg> arg_seg(thread);
	vector<vector<one_index> > seg_index;
	seg_index.resize(thread);
	vector<vector<uint64_t> > read_distance;
	read_distance.resize(thread);
	vector<uint64_t> index_count(thread);
	for(int32_t i=0;i<thread;i++){
		sem_wait(&block_sem);
		usleep(5000);
		
		arg_seg.at(i).l_mem=total_pair.l_mem + mark;
		if(mark + split_count > total_count){
			arg_seg.at(i).count=total_count - mark;
			mark=total_count;
		}else{
			int x=0;
			while(1){
				uint64_t old_bar=total_pair.l_mem[mark+split_count+x-1].bar;//last bar 
				uint64_t new_bar=total_pair.l_mem[mark+split_count+x].bar;//next bar
				if(old_bar != new_bar){
					break;
				}
				x++;
			}
			
			arg_seg.at(i).count=split_count+x;
			mark += (split_count+x);
		}
		
		arg_seg.at(i).outdir=dir;
		arg_seg.at(i).p_index=&(seg_index.at(i));
		arg_seg.at(i).r_index=&(read_distance.at(i));
		arg_seg.at(i).bar_th=bar_th;
		arg_seg.at(i).seg_th=seg_th;
		arg_seg.at(i).gap_th=gap_th;
		arg_seg.at(i).header=header;
		arg_seg.at(i).th_index=i;
		arg_seg.at(i).f_offset=index_count.data()+i;
		
		int res=pthread_create(&tid.at(i),&attr,seg_process,(void*)(arg_seg.data()+i));
		usleep(5000);
		if(res!=0){
			cerr<<"pthread_create error:error code="<<res<<endl;
			exit(1);
		}
	}
	
	while(1){
		int value;
		sem_getvalue(&block_sem,&value);
		if(value == thread){
			break;
		}
	}
	
	bam_hdr_destroy(header);
	cout<<"Merging..."<<endl;
	string cmd;
	string outfile=dir + "/" + prefix + ".sbf";
	cmd="rm -rf "+outfile;
	if(run(cmd.c_str())){
		exit(1);
	}
	
	for(int32_t i=0;i<thread;i++){
		string tempfile=dir+"/segment_"+boost::lexical_cast<string>(i)+".tmp";
		cmd="cat " + tempfile + " >> "+ outfile;
		if(run(cmd.c_str())){
			exit(1);
		}
		
		cmd="rm -rf "+tempfile;
		if(run(cmd.c_str())){
			exit(1);
		}
	}
	
	cout<<"Indexing..."<<endl;
	string outindex=outfile+".bfi";
	cmd="rm -rf "+outindex;
	if(run(cmd.c_str())){
		exit(1);
	}
	ofstream iout(outindex.c_str(),ios::binary);
	uint32_t line_len=sizeof(uint32_t)*4+sizeof(uint64_t)+32;
	uint64_t init_off=0;
	uint64_t switch_bar=0;
	for(int32_t i=0;i<thread;i++){
		for(uint64_t j=0;j < seg_index.at(i).size();j++){
			if(j == 0 && seg_index.at(i).at(j).bar == switch_bar){
				total_seg_bar--;
				continue;
			}
			uint64_t real_off=init_off+seg_index.at(i).at(j).offset *line_len;
			iout.write((char*) &(seg_index.at(i).at(j).bar),sizeof(uint64_t));
			iout.write((char*) &real_off,sizeof(uint64_t));
		}
		init_off +=index_count.at(i)*line_len;
		if(seg_index.at(i).size() > 0){
			switch_bar=seg_index.at(i).back().bar;
		}
	}
	iout.close();
	
	cout<<"Get the LFR read segments...Done"<<endl;
	cout<<"Output statistical info..."<<endl;
	string statfile=dir + "/" + prefix+".stat";
	string gapfile=dir + "/" + prefix+".gap";
	
	ofstream stat(statfile.c_str(),ios::out);
	stat << total_valid_read <<"\t" << total_seg_read << "\t" << total_seg_bar <<endl;
	stat.close();
	
	ofstream gap(gapfile.c_str(),ios::out);
	for(int32_t i=0;i<thread;i++){
		for(uint64_t j=0;j < read_distance.at(i).size();j++){
			gap<<read_distance.at(i).at(j)<<endl;
		}
	}
	gap.close();
	
	cout<<"Output statistical info...Done"<<endl;
	m_end=time(NULL);
	cout<< "Elapsed time " << m_end - m_start << " seconds" << endl;
	total_time=m_end-total_time;
	cout<<"Total elapsed time "<<setprecision(2)<<(double)total_time/3600 << " hours" << endl;
}

void* split_process(void* arg)
{
	sp_arg* sin_arg=(sp_arg*) arg;
	string infile=sin_arg->bam_file;
	string index_file=infile+ ".bai";
	char filt_flag =sin_arg->flag;
	uint32_t is_th=sin_arg->is_th;
	uint64_t valid_read=0;
	
	
	bam_hdr_t *header=NULL;
	hts_itr_t *iter=NULL;
	hts_idx_t *idx=NULL;
	samFile *in=NULL;
	// samFile *out=NULL;
	bam1_t *read = bam_init1();
	in = sam_open(infile.c_str(), "r");
	// out= sam_open(outfile.c_str(), "w");
	header = sam_hdr_read(in);
	idx=sam_index_load(in,index_file.c_str());
	
	for(uint32_t i=0;i<(sin_arg->chrom)->size();i++){
		string region=(header->target_name)[(sin_arg->chrom)->at(i).index_name];
		iter  = sam_itr_querys(idx, header, region.c_str());
		while ( sam_itr_next(in, iter, read) >= 0){
			vector<string> tempname;
			vector<string> tempbar;
			string name=(char * ) read->data;
			// short cnul=(short) read->core.l_extranul;
			// short qlen=(short) read->core.l_qname;
			// char newname[qlen];
			
			boost::split(tempname, name, boost::is_any_of("#"), boost::token_compress_off);
			if(tempname.size() != 2 ){
				cerr <<"Wrong LFR read name format!!"<<endl;
				exit(1);
			}
			if(tempname.at(1) == "0_0_0" ){
				continue;
			}
			boost::split(tempbar, tempname.at(1), boost::is_any_of("_"), boost::token_compress_off);
			name=tempbar.at(0)+"\t"+tempbar.at(1)+"\t"+tempbar.at(2)+"\t"+tempname.at(0);
			
			uint32_t chra=read->core.tid;
			uint32_t chrb=read->core.mtid;
			int32_t posa=read->core.pos;
			int32_t posb=read->core.mpos;
			uint32_t is=abs(read->core.isize);
			uint16_t flag=read->core.flag;
			
			
			//only mapped pair and read1
			if((flag& 0x4) || (flag& 0x8) || (flag& 0x400) || (flag& 0x80)){
				continue;
			}
			
			if(chra != chrb){
				continue;
			}
			
			// if(chra >= 0x100){
				// continue;
			// }
			
			if(is >= 3*is_th){
				continue;
			}
			
			//filter the chromesome name for human
			string pattern="^(?:chr)?[1234567890XYM]{1,2}$";
			regex r (pattern);
			if(filt_flag == 'Y'){
				string tempname=(header->target_name)[chra];
				if(!(regex_match(tempname,r))){
					continue;	
				}
			}
			
			read_pair onepair;
			onepair.bar=covertbar(tempbar);
			onepair.index_name= chra;
			onepair.mid_pos=(uint32_t)(posa+1+posb+1)/2;
			valid_read++;
			
			//check the size
			if((((sin_arg->mem)->size)-((sin_arg->mem)->count)) < 100){
				(sin_arg->mem)->l_mem=(read_pair *) realloc((sin_arg->mem)->l_mem,((uint64_t)((sin_arg->mem)->size * 1.5)) *  sizeof(read_pair) );
				(sin_arg->mem)->size=(uint64_t)((sin_arg->mem)->size * 1.5);
			}
			
			((sin_arg->mem)->l_mem)[(sin_arg->mem)->count]=onepair;
			((sin_arg->mem)->count)++;
		}	
	}
	sam_close(in);

	bam_hdr_destroy(header);
	hts_itr_destroy(iter);
	hts_idx_destroy(idx);
	bam_destroy1(read);
	
	pthread_mutex_lock(&mutex1);
	(sin_arg->tid_id)->push_back(sin_arg-> th_index);
	total_valid_read+=valid_read;
	pthread_mutex_unlock(&mutex1);
	sem_post(&block_sem);
	return NULL;
}

void* seg_process(void* arg)
{
	seg_arg* sin_arg=(seg_arg*) arg;
	uint32_t prefix=sin_arg->th_index;
	read_pair * p=sin_arg->l_mem;
	vector<one_index> * p_index=sin_arg->p_index;
	vector<uint64_t> * r_index=sin_arg->r_index;
	uint64_t count=sin_arg->count;
	string dir=sin_arg->outdir;
	uint32_t bar_th=sin_arg->bar_th;
	uint32_t seg_th=sin_arg->seg_th;
	uint32_t gap_th=sin_arg->gap_th;
	bam_hdr_t *header=sin_arg->header;
	uint64_t * final_off=sin_arg->f_offset;
	uint64_t seg_bar=0;
	uint64_t seg_read=0;
	
	string outfile=dir+"/segment_"+boost::lexical_cast<string>(prefix)+".tmp";
	ofstream fout(outfile.c_str(),ios::binary);
	vector<read_pair> bar;
	uint64_t last_bar=0;
	one_index last_index;
	last_index.bar=0;
	last_index.offset=0;
	for(uint64_t i=0;i<count;i++){
		read_pair one=p[i];
		if(one.bar != last_bar){
			if(bar.size() >= bar_th){
				vector<seg> segment;
				onebar_process(bar,segment,seg_th,gap_th,seg_read,r_index);
				if(segment.size() > 0){
					seg_bar++;
					//check the bar,update the index
					uint64_t tempbar=segment.front().bar;
					tempbar= tempbar >> 20;
					if(last_index.bar != tempbar){
						last_index.bar=tempbar;
						p_index->push_back(last_index);
					}
					//output the segment
					for(uint32_t j=0;j< segment.size();j++){
						char name[32];
						memcpy(name,(header->target_name)[segment.at(j).index_name],31*sizeof(char));
						name[31]='\0';
						fout.write((char*) &(segment.at(j).bar),sizeof(uint64_t));
						fout.write((char*) &(segment.at(j).seg_index),sizeof(uint32_t));
						fout.write(name,32);
						fout.write((char*) &(segment.at(j).s_pos),sizeof(uint32_t));
						fout.write((char*) &(segment.at(j).e_pos),sizeof(uint32_t));
						fout.write((char*) &(segment.at(j).sup_count),sizeof(uint32_t));
						last_index.offset++;
					}
				}
			}
			bar.clear();
			last_bar=one.bar;
		}
		bar.push_back(one);	
	}
	
	if(bar.size() >= bar_th){
		vector<seg> segment;
		onebar_process(bar,segment,seg_th,gap_th,seg_read,r_index);
		if(segment.size() > 0){
			seg_bar++;
			//check the bar,update the index
			uint64_t tempbar=segment.front().bar;
			tempbar= tempbar >> 20;
			if(last_index.bar != tempbar){
				last_index.bar=tempbar;
				p_index->push_back(last_index);
			}
			//output the segment
			for(uint32_t j=0;j< segment.size();j++){
				char name[32];
				memcpy(name,(header->target_name)[segment.at(j).index_name],31*sizeof(char));
				name[31]='\0';
				fout.write((char*) &(segment.at(j).bar),sizeof(uint64_t));
				fout.write((char*) &(segment.at(j).seg_index),sizeof(uint32_t));
				fout.write(name,32);
				fout.write((char*) &(segment.at(j).s_pos),sizeof(uint32_t));
				fout.write((char*) &(segment.at(j).e_pos),sizeof(uint32_t));
				fout.write((char*) &(segment.at(j).sup_count),sizeof(uint32_t));
				last_index.offset++;
			}
		}
	}
	fout.close();
	*final_off=last_index.offset;
	pthread_mutex_lock(&mutex1);
	total_seg_bar+=seg_bar;
	total_seg_read+=seg_read;
	pthread_mutex_unlock(&mutex1);
	sem_post(&block_sem);
	return NULL;
}

void onebar_process(vector<read_pair> &bar, vector<seg> &segment, uint32_t seg_th, uint32_t gap_th, uint64_t &seg_read, vector<uint64_t> * &r_index)
{
	uint64_t barcode=bar.at(0).bar;
	vector<read_pair> one_seg;
	int64_t last_name=-1;
	uint32_t last_pos=0;
	uint32_t seg_count=0;
	seg one_segment;
	for(uint32_t i=0; i< bar.size();i++){
		if((int64_t) bar.at(i).index_name != last_name){
			if(one_seg.size() >= seg_th){
				seg_read+=one_seg.size();
				one_segment.bar=barcode;
				one_segment.seg_index=seg_count;
				one_segment.index_name=one_seg.front().index_name;
				one_segment.s_pos=one_seg.front().mid_pos;
				one_segment.e_pos=one_seg.back().mid_pos;
				one_segment.sup_count=one_seg.size();
				segment.push_back(one_segment);
				seg_count++;
			}
			one_seg.clear();
			
			one_seg.push_back(bar.at(i));
			last_name=bar.at(i).index_name;
			last_pos=bar.at(i).mid_pos;
		}else{
			uint64_t dis=bar.at(i).mid_pos - last_pos;
			r_index->push_back(dis);
			if(dis >= gap_th){
				if(one_seg.size() >= seg_th){
					seg_read+=one_seg.size();
					one_segment.bar=barcode;
					one_segment.seg_index=seg_count;
					one_segment.index_name=one_seg.front().index_name;
					one_segment.s_pos=one_seg.front().mid_pos;
					one_segment.e_pos=one_seg.back().mid_pos;
					one_segment.sup_count=one_seg.size();
					segment.push_back(one_segment);
					seg_count++;
				}
				one_seg.clear();
			}
			last_pos=bar.at(i).mid_pos;
			one_seg.push_back(bar.at(i));	
		}
	}
	
	if(one_seg.size() >= seg_th){
		seg_read+=one_seg.size();
		one_segment.bar=barcode;
		one_segment.seg_index=seg_count;
		one_segment.index_name=one_seg.front().index_name;
		one_segment.s_pos=one_seg.front().mid_pos;
		one_segment.e_pos=one_seg.back().mid_pos;
		one_segment.sup_count=one_seg.size();
		segment.push_back(one_segment);
		seg_count++;
	}
	return;
}

bool comp_chrom(const chr &a,const chr &b)
{
	return a.tid_len > b.tid_len;
}

uint64_t covertbar(vector<string> & tempbar)
{
	uint64_t i=0;
	i=boost::lexical_cast<uint64_t>(tempbar.at(0));
	i=i<<20;
	i|=boost::lexical_cast<uint64_t>(tempbar.at(1));
	i=i<<20;
	i|=boost::lexical_cast<uint64_t>(tempbar.at(2));
	return i;
}

bool comp_read(const read_pair &a,const read_pair &b)
{
	if(a.bar != b.bar){
		return a.bar < b.bar;	
	}else if(a.index_name != b.index_name ){
		return a.index_name < b.index_name;
	}else{
		return a.mid_pos < b.mid_pos;	
	}
}

int run(const char *command)
{
    int ret = system(command);
    if (-1 == ret)
    {
        fprintf(stderr,"sh_shell fork fail");
        return -1;
    }
    else  
    {  
        if (WIFEXITED(ret))  
        {
            if (0 == WEXITSTATUS(ret))  
            {   
                return 0;
            }  
            else  
            {
                fprintf(stderr,"sh_shell fail, shell status: %d", WEXITSTATUS(ret));
                return -1;
            }  
        }  
        else  
        {  
            fprintf(stderr,"sh_shell exit status = %d", WEXITSTATUS(ret));  
            return -1;
        }  
    } 
}