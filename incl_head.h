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

	const string word="what starts here changes the world";	
	const int word_length=4;
	const int L=4;
	const int KG=64,Kopt=3;  // from 00000(means space) 00001(a) ... 11010(z)
	const int J=2*KG*L+1,Tseq=word_length*L;
	const int update_num=4000, Inner_num=30;
	const double mu=0.3,cost_un=1.0,cost_mis=1.5;
	const double prefer=0.2;