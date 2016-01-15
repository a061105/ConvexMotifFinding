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

	const string word="aaaaaaaaacccccdddp";	
	const int word_length=word.size();
	const int L=4;
	const int KG=4,Kopt=4;
	const int J=2*KG*L+1,Tseq=word_length*L;
	const int update_num=4000, Inner_num=30;
	const double mu=0.3,cost_un=1.0,cost_mis=1.5;