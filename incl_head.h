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

	const string word="veni vidi vici";
	const int word_length=14,space_num=1;
	const int L=5,Lmin=4;
	const int KG=32,Kopt=7;  // from 00000(means space) 00001(a) ... 11010(z)
	const int J=2*KG*L+1,Tseq=word_length*Lmin;
	const int update_num=4000, Inner_num=30;
	const double mu=0.3,cost_un=1.0,cost_mis=1.5;
	const double prefer=0.1,global_prefer=prefer/KG/L;
	const double inner_eps=1e-5, outer_eps=1e-5;