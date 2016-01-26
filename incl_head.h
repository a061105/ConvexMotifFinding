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

	const string word="aaaa";
	const int word_length=4,space_num=1;
	const int L=5,Lmin=4;
	const int KG1=32,KG2=16,Kopt=2;  // from 00000(means space) 00001(a) ... 11010(z)
	const int J1=32*2*L+1,J2=16*2*Lmin,Tseq=word_length*Lmin+space_num;
	const int update_num=4000, Inner_num=30;
	const double mu=0.3,cost_un=1.0,cost_mis=1.5;
	const double prefer=0.1,short_prefer=0*prefer/(KG1+KG2);
	const double inner_eps=1e-5, outer_eps=1e-5;