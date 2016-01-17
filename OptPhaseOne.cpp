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
extern vector<double> OptStep2( MAT_D W1,
						MAT_D W2,
						MAT_D YY,
						MAT_D DIR,
						vector<int> pattern_ind);
extern vector<int> Frank_Wolfe(	MAT_D Gradf);
extern void Report_dir(vector<int> AT);
extern MAT_D GroupConsDir(MAT_D Gradw2, vector<int>& pattern_ind);

//----------------------Update W1------------------------------
MAT_D OptPhaseOne(	MAT_D CC,
					MAT_D WW1,
					MAT_D WW2,
					MAT_D YY,
					int num_k)
{
	
	MAT_D Gradw1;
	vector<int> Aton;
	double Step_size=2.0/(num_k+2);
	
	Gradw1=Gradf(CC,WW1,WW2,YY);
	Aton=Frank_Wolfe(Gradw1);
	Step_size=OptStep1(CC,WW1,WW2,YY,Aton);
	if(Step_size>1) Step_size=1.0;
	if(Step_size<0) Step_size=0;
	//cout<<"Step_size: "<<Step_size<<endl;
	//Report_dir(Aton);
	// update W1
	for(int t=0; t<Tseq; t++){
		for(int j=0; j<J; j++){
			WW1[j][t]*=(1.0-Step_size);
		}
	}
	for(int t=0; t<Tseq; t++){
		int j=Aton[t];
		WW1[j][t]+=Step_size;
	}
	return WW1;
}

//----------------------Update W2------------------------------
MAT_D OptPhaseTwo(	MAT_D CC,
					MAT_D WW1,
					MAT_D WW2,
					MAT_D YY,
					int num_k)
{
	
	MAT_D Gradw2;
	MAT_D W2dir;
	vector<int> pattern_ind(Kopt,0);
	vector<double> Step_size(KG,2.0/(num_k+2));
	Gradw2=Gradf(WW1,WW2,YY);
	W2dir=GroupConsDir(Gradw2, pattern_ind);
	Step_size=OptStep2(WW1,WW2,YY,W2dir,pattern_ind);
	//cout<<"Step_size: "<<Step_size<<endl;
	//Report_dir(Aton);
	// update W2
	for(int t=0; t<Tseq; t++){
		for(int pat=0; pat<Kopt; pat++){
			int kk=pattern_ind[pat];
			for(int j=2*L*kk+1; j<2*L*(kk+1)+1; j++){
			WW2[j][t]=WW2[j][t]*(1.0-Step_size[kk])+W2dir[j][t]*Step_size[kk];
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
	return WW2;
}

//Compute Original Grad for W1
MAT_D Gradf(	MAT_D CC,
				MAT_D WW1,
				MAT_D WW2,
				MAT_D YY)
{
	MAT_D GG=CC;
	for(int i=0; i<Tseq; i++){
		for(int j=0; j<J; j++){
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
	MAT_D GG(J,col_dou);
	for(int i=0; i<Tseq; i++){
		for(int j=0; j<J; j++){
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
	for(int t=0; t<Tseq; t++){
		for(int j=0; j<J; j++){
			YR[j][t]=YR[j][t]*(1.0)+mu*(WW1[j][t]-WW2[j][t]);
		}
	}
	return YR;
}

