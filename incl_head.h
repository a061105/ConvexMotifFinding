#ifndef FRACTION_H
#define FRACTION_H
#endif
#include<iostream>
#include<time.h>
#include<cmath>
#include<string>
#include<fstream>
#include<vector>
#include<stdio.h>
#include<stdlib.h>
#include<algorithm>
#include <functional>
#define random  (rand()/(double)(RAND_MAX)); 
#define MAX_NUMBER 999999
using namespace std;
typedef vector<vector<double> > MAT_D; 

	//const string word="aaab";	
	const string InSeq="000010000100001000011";
	const int word_length=4;
	const int L=6,Lmin=5;
	const int KG=32,Kopt=2;  // from 00000(means space) 00001(a) ... 11010(z)
	const int J=2*KG*L+1,Tseq=InSeq.size();
	const int update_num=4000, Inner_num=30;
	const double mu=0.3,cost_un=1.0,cost_mis=1.5;
	const double prefer=0.1,global_prefer=prefer/KG;
	const double penalty_imcomplete=0.0;
	const double inner_eps=0.01, outer_eps=1e-5;