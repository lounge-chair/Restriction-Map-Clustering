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


int overlap(int s1,int e1, int s2, int e2){

    if(s1<=s2 && e1>=e2){
        return (e2-s2);
    }
    else if(s2<=s1 && e2>=e1){
        return (e1-s1);
    }
    else if(s1<=s2 && e1<=e2){
        return (e1-s2);
    }
    else if(s2<=s1 && e1>=e2){
        return (e2-s1);
    }
    else{
        cout<<"weird case"<<endl;
        return 0;
    }


}

int main(){

    ifstream infofile("all_info");

    vector<int> orientation;
    vector<int> starts,ends;
    vector<int> cluster_assignment;
    vector<string> rmap_belong;


    map<string,int> nameid;
    int lastname=0,rmap_cnt=0;

    ifstream valfile("ecoli_cov300.val");

    string str;

    while(getline(valfile,str)){
        stringstream ss11(str);

        ss11>>str;
        if(nameid.find(str)==nameid.end()){
            nameid[str]=lastname++;
        }

        getline(valfile,str);

        stringstream ss12(str);

        string dummy;
        float f1;

        ss12>>dummy>>dummy;

        ss12>>f1;

        rmap_cnt++;

        getline(valfile,str);
    }

    string sttt;

    while(getline(infofile,sttt)){

        stringstream sstt(sttt);

        string dum;
        int len;
        sstt>>dum>>len>>dum;

        int dumin;

        sstt>>dumin;
        int orr;
        sstt>>dum>>orr;
        ends.push_back(dumin+len);
        sstt>>orr;

        starts.push_back(dumin);
        orientation.push_back(orr);

    }


    ifstream clusterfile("/ufrc/boucher/kingdgp/clustering/myclustering_k6_bin1600/model_predictions.csv");

    int maxcluster=0;

    while(getline(clusterfile,sttt)){

        std::size_t found = sttt.find(",");

        sttt=sttt.substr(found+1);

        stringstream stss(sttt);

//        cout<<sttt<<endl;
        int cl;

        stss>>cl;
        cluster_assignment.push_back(cl);

        if(cl>maxcluster)
            maxcluster=cl;

    }

    vector < vector <int> > rmap_to_cluster;

    vector < vector <int> > cluster_to_rmap;

    vector<int> temp;

    cluster_to_rmap.resize(maxcluster,temp);

    ifstream rmapfile("6merinfo.txt");

    int lineno=0;

    rmap_to_cluster.resize(starts.size(),temp);

    while(getline(rmapfile,sttt)){

        rmap_belong.push_back(sttt);

        rmap_to_cluster[nameid.find(sttt)->second].push_back(cluster_assignment[lineno]);

        lineno++;

    }


    for(int i=0;i<cluster_assignment.size();i++){

        cluster_to_rmap[cluster_assignment[i]-1].push_back(nameid.find(rmap_belong[i])->second);


    }


    ofstream cltorm("cluster_rmap.txt");
    ofstream rmtocl("rmap_to_cl.txt");

    for(int i=0;i<rmap_to_cluster.size();i++){

        rmtocl<<i<<" : ";

        for(int j=0;j<rmap_to_cluster[i].size();j++){

            rmtocl<<rmap_to_cluster[i][j]<<" ";
        }
        rmtocl<<endl;
    }

    for(int i=0;i<cluster_to_rmap.size();i++){

        cltorm<<i<<" : ";

        for(int j=0;j<cluster_to_rmap[i].size();j++){

            cltorm<<cluster_to_rmap[i][j]<<" ";
        }
        cltorm<<endl;
    }


    vector< vector<int> > relation_matrix;


    vector<int> temvec(rmap_to_cluster.size(),0);
    relation_matrix.resize(rmap_to_cluster.size(),temvec);

    vector< vector <int> > precision_true_relations;
    vector< vector <int> > recall_true_relations;

    precision_true_relations.resize(rmap_to_cluster.size(),temp);
    recall_true_relations.resize(rmap_to_cluster.size(),temp);

    for(int i=0;i<rmap_to_cluster.size();i++){

        if(rmap_to_cluster[i].size()==0){
                continue;
        }

        for(int j=0;j<rmap_to_cluster.size();j++){
            if(j==i || rmap_to_cluster[j].size()==0 || orientation[i]!=orientation[j])
                continue;

            if(overlap(starts[i],ends[i],starts[j],ends[j])>10000)
                precision_true_relations[i].push_back(j);

            if(overlap(starts[i],ends[i],starts[j],ends[j])>70000)
                recall_true_relations[i].push_back(j);



        }

    }

    ofstream pretrurel("precision_true_relations.txt");
    ofstream rectrurel("recall_true_relations.txt");




    ofstream rmaprels("Rmap_relations.txt");

    cout<<endl<<relation_matrix.size()<<" x "<<relation_matrix[0].size()<<endl;

    vector< vector <int> > found_rels(relation_matrix.size(),temp);


    for(int i=0;i<rmap_to_cluster.size();i++){

        for(int j=0;j<rmap_to_cluster[i].size();j++){

            int cluster=rmap_to_cluster[i][j]-1;


            for(int k=0;k<cluster_to_rmap[cluster].size();k++){

                relation_matrix[i][cluster_to_rmap[cluster][k]]++;
            }


        }
        rmaprels<<i<<" : ";
        for(int j=0;j<relation_matrix[i].size();j++){

            if(relation_matrix[i][j]>10){
                found_rels[i].push_back(j);
                rmaprels<<j<<","<<relation_matrix[i][j]<<"\t";
            }


        }

        rmaprels<<endl;


    }


    int correct=0,incorrect=0,totcorrect=0,totfalse=0,recor=0,refal=0;

    ofstream prefile("precision_file.txt");

    for(int i=0;i<precision_true_relations.size();i++){
        pretrurel<<i<<" : ";
        correct=0;incorrect=0;

        vector<int> track(precision_true_relations.size(),0);

        for(int j=0;j<precision_true_relations[i].size();j++){
            pretrurel<<precision_true_relations[i][j]<<" ";
            track[precision_true_relations[i][j]]++;
        }

        for(int j=0;j<found_rels[i].size();j++){
            if(track[found_rels[i][j]]>0)
                correct++;
            else
                incorrect++;
        }

        totcorrect+=correct;
        totfalse+=incorrect;

        prefile<<i<<" : "<<correct<<" "<<incorrect<<" "<<found_rels[i].size()<<endl;

        pretrurel<<endl;

        rectrurel<<i<<" : ";

        track.resize(precision_true_relations.size(),0);

        for(int j=0;j<recall_true_relations[i].size();j++){
            rectrurel<<recall_true_relations[i][j]<<" ";
            track[recall_true_relations[i][j]]++;
        }
        rectrurel<<endl;

        correct=0;incorrect=0;

        for(int j=0;j<found_rels[i].size();j++){
            if(track[found_rels[i][j]]>0)
                correct++;

        }
        incorrect=recall_true_relations[i].size()-correct;

        recor+=correct;
        refal+=incorrect;

    }


    cout<<"Precision: "<<(float)totcorrect*100/(totcorrect+totfalse)<<endl;
    cout<<"Recall: "<<(float)recor*100/(recor+refal)<<endl;



}
