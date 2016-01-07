#ifndef FRACTION_H
#define FRACTION_H
#endif
#include<iostream>
#include<cmath>
#include<string>
#include<vector>
#include<stdio.h>
#include<stdlib.h>
#define random ((rand()%100)/100.0)
#define MAX_NUMBER 999999
using namespace std;
typedef vector<vector<double> > MAT_D; 

	const int word_length=8,L=8;
	const int KG=3,Kopt=3;
	const int J=2*KG*L+1,Tseq=word_length*L;
	const int update_num=100, Inner_num=30;
	const double mu=0.1,cost_un=1.0,cost_mis=1.0;