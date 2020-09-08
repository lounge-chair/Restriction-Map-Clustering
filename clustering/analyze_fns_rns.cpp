#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <map>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm>
#include <set>
#include <ctime>
#include <vector>
#include <stdlib.h>
#include <typeinfo>
#include <string>

#define NDEBUG
#include <cassert>

#define Pi 3.14159265
#define pi 3.14156

// d as genomic distance

using namespace std;

float threshold_rep_nodes=40;
int rangefut=8;
int multiplicity=5;

vector <int> covered;

vector<int> fun_vec_read;

vector< vector<float> > real_values;

vector<int> fun_vec_write;

vector <string> neighs;

vector< set < int > > path_tracker;

vector< vector <int> > final_contigs;

vector< vector < int > > future_of_contigs;

vector < vector< set < int > > > complete_future;


vector< vector <float> > final_contigs_real;

vector< vector <int> > final_contig_omgrams;

vector <int>  lng_contig;
vector <float>  lng_contig_real;

vector< int > long_contig_kmers;


set < int > futures;

set < int > start_nodes;

vector <int> nodes_in_path;

vector< vector<float> > not_sure;


vector<float> this_stitch;
vector<int> this_stitch_nodes;

vector< vector<int> > contig_nodes;
vector<bool> nodes_covered_start;
vector<bool> nodes_covered_ctnd;

//vector < pair<int ,int> > optimized_overlap_alignment(vector< float >& rmap1, vector< float >& rmap2);

int min(int a,int b){

    if(a<=b)
        return a;
    else
        return b;
}

int max(int a,int b){

    if(a>=b)
        return a;
    else
        return b;
}


int main(){

    ifstream infile("fns_rns.rf");

    ofstream rmaout("reference.out");

    string str;

    getline(infile,str);

    getline(infile,str);

    stringstream sss(str);

    vector<int> fwd,rev,locf,locr;
    int tem;

    while(sss>>tem){
        locf.push_back(tem);
    }

    getline(infile,str);

    getline(infile,str);


    stringstream sss2(str);

    while(sss2>>tem){
        locr.push_back(tem);
    }

    for(int i=1;i<locf.size();i++){
        tem=locf[i]-locf[i-1];

        fwd.push_back(tem);

    }

    for(int i=1;i<locr.size();i++){
        tem=locr[i]-locr[i-1];

        rev.push_back(tem);

    }


    ofstream outfwd("analysisout.fwd");

    outfwd<<locf[0]<<endl;

    for(int i=0;i<fwd.size();i++){

        outfwd<<"\t"<<fwd[i]<<endl;

        outfwd<<locf[i+1]<<endl;
    }

    ofstream outrev("analysisout.rev");

    outrev<<locr[0]<<endl;

    for(int i=0;i<rev.size();i++){

        outrev<<"\t"<<rev[i]<<endl;

        outrev<<locr[i+1]<<endl;
    }

    rmaout<<"forward_reference"<<endl<<" \t BspQI \t BspQI \t";
    for(int i=0;i<fwd.size();i++){
        rmaout<<(float)fwd[i]/1000<<" ";

    }

    rmaout<<endl<<endl;

    rmaout<<"reverse_reference"<<endl<<" \t BspQI \t BspQI \t";
    for(int i=0;i<rev.size();i++){
        rmaout<<(float)rev[i]/1000<<" ";

    }

    rmaout<<endl<<endl;



}
