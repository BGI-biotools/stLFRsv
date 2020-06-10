//sort program for LFR data 
//by guojunfu
//lascqy@hotmail.com
#include <string>
#include <pthread.h>
#include <vector>
#define INIT_SPLIT 30000000
using namespace std;
//for multi thread;
sem_t block_sem;
pthread_mutex_t mutex1;
// pthread_mutex_t mutex2;
uint64_t total_valid_read=0;
uint64_t total_seg_read=0;
uint64_t total_seg_bar=0;

struct chr
{
	uint32_t tid_len;
	uint32_t index_name;
};

struct read_pair
{
	uint64_t bar;
	uint32_t mid_pos;
	uint32_t index_name;
};

struct mem_block
{
	uint64_t size;
	uint64_t count;
	read_pair * l_mem;
};

struct one_index
{
	uint64_t bar;
	uint64_t offset;
};

struct seg
{
	uint64_t bar;
	uint32_t seg_index;
	uint32_t s_pos;
	uint32_t e_pos;
	uint32_t sup_count;
	uint32_t index_name;
};

struct sp_arg
{
	int th_index;
	vector<chr>* chrom;
	mem_block * mem;
	vector<int>* tid_id;
	char flag;
	uint32_t is_th;
	string bam_file;
};

struct seg_arg
{
	uint64_t count;
	uint32_t th_index;
	uint32_t bar_th;
	uint32_t seg_th;
	uint32_t gap_th;
	uint64_t * f_offset;
	bam_hdr_t * header;
	read_pair * l_mem;
	vector<one_index> * p_index;
	vector<uint64_t> * r_index;
	string outdir;
};

void* split_process(void* arg);
void* seg_process(void* arg);
void onebar_process(vector<read_pair> &bar, vector<seg> &segment, uint32_t seg_th, uint32_t gap_th, uint64_t &seg_read, vector<uint64_t> * &r_index);
bool comp_chrom(const chr &a,const chr &b);
bool comp_read(const read_pair &a,const read_pair &b);
uint64_t covertbar(vector<string> & tempbar);
int run(const char *command);