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
#include <bits/stdc++.h>
#include <climits>



#define NDEBUG
#include <cassert>

#define Pi 3.14159265
#define pi 3.14156

using namespace std;

int main(int argc, char *argv[]){

    if(argc<3){
	cout<<"usage: ./cometrel <kmer_file> <prediction_outputfile>"<<endl;
	return(1);}

    string str;
    ifstream inkmer(argv[1]);
    ofstream outfile(argv[2]);
    outfile<<"kmer_index,bucket_index"<<endl;
    float    bucket=4;

    vector<string> kmers;
    map<string,int> mapkmer;
    int lastmap=1;
    int kmercnt=1;

//    ofstream check("check");

    while (getline(inkmer,str)){

        istringstream ss2(str);
        string thiskmer;
        float input;
        stringstream outst;

        while(ss2>>input){
            float infl= input/bucket;
            int temin = round(infl);
//            check<<input<<"/"<<infl<<"/"<<temin<<" ";
            outst<<temin<<"_";
        }
//        check<<endl;
        thiskmer = outst.str();

        if(mapkmer.find(thiskmer)!=mapkmer.end()){
            outfile<<kmercnt++<<"/"<<thiskmer<<","<<mapkmer.find(thiskmer)->second<<endl;
        }
        else{
            mapkmer[thiskmer]=lastmap++;
            outfile<<kmercnt++<<"/"<<thiskmer<<","<<mapkmer.find(thiskmer)->second<<endl;
        }

    }

    cout<<"Total number of kmers: "<<kmercnt<<endl;

    cout<<"Total number of kmer_buckets: "<<lastmap<<endl;
}
