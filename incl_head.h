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
#include <fstream>
#define random  (rand()/(double)(RAND_MAX)); 
#define MAX_NUMBER 999999
using namespace std; 
typedef vector<vector<double> > MAT_D; 

	const string word="veni vidi vici";
	const int word_length=14,space_num=2;
	const int L=7,Lmin=4;
	const int KG1=128,KG2=16,Kopt=7;  // from 00000(means space) 00001(a) ... 11010(z)
	const int J1=KG1*2*L+1,J2=KG2*2*Lmin,Tseq=word_length*L-3*space_num;
	const int update_num=1000, Inner_num=30;
	const double mu=0.3,cost_un=1.0,cost_mis=1.5;
	const double prefer=0.1,short_prefer=prefer/(L+Lmin);
	const double inner_eps=1e-5, outer_eps=1e-5;