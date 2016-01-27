#include"incl_head.h"
MAT_D Gradf(		MAT_D CC,
					MAT_D WW1,
					MAT_D WW2,
					MAT_D YY);
MAT_D Gradf(		MAT_D WW1,
					MAT_D WW2,
					MAT_D YY);
extern double OptStep1( MAT_D CC,
						MAT_D W1,
						MAT_D W2,
						MAT_D YY,
						vector<int> Aton);
extern double OptStep2( MAT_D W1,
						MAT_D W2,
						MAT_D YY,
						MAT_D DIR,
						vector<int> pattern_ind);
extern vector<int> Frank_Wolfe(	MAT_D Gradf);
extern void Report_dir(vector<int> AT);
extern MAT_D GroupConsDir(MAT_D Gradw2, vector<int>& pattern_ind);

//----------------------Update W1------------------------------
void OptPhaseOne(	MAT_D& CC,
					MAT_D& WW1,
					MAT_D& WW2,
					MAT_D& YY,
					double& Step_size)
{
	
	MAT_D Gradw1;
	vector<int> Aton;
	Gradw1=Gradf(CC,WW1,WW2,YY);
	Aton=Frank_Wolfe(Gradw1);
	Step_size=OptStep1(CC,WW1,WW2,YY,Aton);
	if(Step_size>1) Step_size=1.0;
	if(Step_size<0) Step_size=0;
	//cout<<"Step_size: "<<Step_size<<endl;
	//Report_dir(Aton);
	// update W1
	for(int j=0; j<J1+J2; j++){
		for(int t=0; t<Tseq; t++){
			WW1[j][t]*=(1.0-Step_size);
		}
	}
	for(int t=0; t<Tseq; t++){
		int j=Aton[t];
		WW1[j][t]+=Step_size;
	}
	//cout<<Step_size<<endl;
	return;
}

//----------------------Update W2------------------------------
void OptPhaseTwo(	MAT_D& CC,
					MAT_D& WW1,
					MAT_D& WW2,
					MAT_D& YY,
					double& MaxStep_size)
{
	
	MAT_D Gradw2;
	MAT_D W2dir;
	vector<int> pattern_ind(Kopt,0);
	double Step_size=MaxStep_size;
	Gradw2=Gradf(WW1,WW2,YY);
	W2dir=GroupConsDir(Gradw2, pattern_ind);
	Step_size=OptStep2(WW1,WW2,YY,W2dir,pattern_ind);
	//find max step_size
	MaxStep_size=Step_size;
	//cout<<Step_size<<endl;
	//cout<<endl;
	//Report_dir(Aton);
	// update W2
	
	for(int j=0; j<J1+J2; j++){
			for(int t=0; t<Tseq; t++){
				WW2[j][t]*=(1.0-Step_size);
			}
	}

	
	for(int pat=0; pat<Kopt; pat++){
		int kk=pattern_ind[pat];
		int the_length=L;
		if(kk>=KG1) the_length=Lmin;
		for(int inn=0; inn<2*the_length; inn++){
			int j=2*L*kk+inn+1;
			if(kk>=KG1) j=2*Lmin*(kk-KG1)+J1+inn;
			for(int t=0; t<Tseq; t++){
				WW2[j][t]+=W2dir[j][t]*Step_size;
			}
		}
	}
	// Fit unassigned entrees
	for(int t=0; t<Tseq; t++){
		double W0t=WW1[0][t]+YY[0][t]/mu;
		if(W0t>1) W0t=1.0;
		if(W0t<0) W0t=0;
		WW2[0][t]=W0t;
	}
	//cout<<MaxStep_size<<endl;
	return;
}

//Compute Original Grad for W1
MAT_D Gradf(	MAT_D CC,
				MAT_D WW1,
				MAT_D WW2,
				MAT_D YY)
{
	MAT_D GG=CC;
	for(int j=0; j<J1+J2; j++){
		for(int i=0; i<Tseq; i++){
			GG[j][i]+=mu*(WW1[j][i]-WW2[j][i]+YY[j][i]/mu);
		}
	}
	return GG;
}
// reload Gradf for W2
MAT_D Gradf(	MAT_D WW1,
				MAT_D WW2,
				MAT_D YY)
{
	vector<double> col_dou(Tseq,0.0);
	MAT_D GG(J1+J2,col_dou);
	for(int j=0; j<J1+J2; j++){
		for(int i=0; i<Tseq; i++){
			GG[j][i]-=mu*(WW1[j][i]-WW2[j][i]+YY[j][i]/mu);
		}
	}
	return GG;
}

// Update Y=Y+mu*(W1-W2)
MAT_D UpdateY(	MAT_D WW1,
				MAT_D WW2,
				MAT_D YY){
	MAT_D YR(YY);
	for(int j=0; j<J1+J2; j++){
		for(int t=0; t<Tseq; t++){
			YR[j][t]=YR[j][t]+mu*(WW1[j][t]-WW2[j][t]);
		}
	}
	return YR;
}

void SamplePhaseTwo(		MAT_D& CC,
							MAT_D& WW1,
							MAT_D& WW2,
							MAT_D& YY,
							double& MaxStep_size,
							vector<int>& Aton,
							double& loss){
	
	MAT_D Fake_CC(CC);
	MAT_D Gradw2;
	MAT_D W2dir;
	vector<int> pattern_ind(Kopt,0);
	double Step_size=MaxStep_size;
	Gradw2=Gradf(WW1,WW2,YY);
	W2dir=GroupConsDir(Gradw2, pattern_ind);
	Step_size=OptStep2(WW1,WW2,YY,W2dir,pattern_ind);
	//find max step_size
	MaxStep_size=Step_size;
	//cout<<Step_size<<endl;
	//cout<<endl;
	// add penalty to unchosen patterns
	for(int j=0; j<J1+J2; j++){
			for(int t=0; t<Tseq; t++){
				Fake_CC[j][t]+=cost_mis;
			}
	}
	for(int indi=0; indi<Kopt; indi++){
		int kk=pattern_ind[indi];
		int pat_length=L;
		if(kk>=KG1) pat_length=Lmin;
		for(int j=0; j<2*pat_length; j++){
			int jj=2*L*kk+1+j;
			if(kk>=KG1) jj=2*Lmin*(kk-KG1)+J1+j;
			for(int t=0; t<Tseq; t++){
				Fake_CC[jj][t]-=cost_mis;
			}
		}
	}

	Aton=Frank_Wolfe(Fake_CC);
	loss=0;
	for(int t=0; t<Tseq; t++){
		int j=Aton[t];
		loss+=Fake_CC[j][t];
	}
	// update W2



	

	for(int j=0; j<J1+J2; j++){
			for(int t=0; t<Tseq; t++){
				WW2[j][t]*=(1.0-Step_size);
			}
	}

	
	for(int pat=0; pat<Kopt; pat++){
		int kk=pattern_ind[pat];
		int the_length=L;
		if(kk>=KG1) the_length=Lmin;
		for(int inn=0; inn<2*the_length; inn++){
			int j=2*L*kk+inn+1;
			if(kk>=KG1) j=2*Lmin*(kk-KG1)+J1+inn;
			for(int t=0; t<Tseq; t++){
				WW2[j][t]+=W2dir[j][t]*Step_size;
			}
		}
	}
	// Fit unassigned entrees
	for(int t=0; t<Tseq; t++){
		double W0t=WW1[0][t]+YY[0][t]/mu;
		if(W0t>1) W0t=1.0;
		if(W0t<0) W0t=0;
		WW2[0][t]=W0t;
	}
	//cout<<MaxStep_size<<endl;
	return;
}