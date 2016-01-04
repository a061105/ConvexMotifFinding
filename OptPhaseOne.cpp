#include"incl_head.h"
MAT_D Gradf(		MAT_D CC,
					MAT_D WW1,
					MAT_D WW2,
					MAT_D YY);
MAT_D Gradf(		MAT_D WW1,
					MAT_D WW2,
					MAT_D YY);

extern vector<int> Frank_Wolfe(	MAT_D Gradf);
extern void Report_dir(vector<int> AT);
extern MAT_D GroupConsDir(MAT_D Gradw2);

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
	double Step_size=2.0/(num_k+2);
	Gradw2=Gradf(WW1,WW2,YY);
	W2dir=GroupConsDir(Gradw2);
	//Report_dir(Aton);
	// update W2
	for(int t=0; t<Tseq; t++){
		for(int j=0; j<J; j++){
			WW2[j][t]+=Step_size*W2dir[j][t];
		}
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
			GG[j][i]+=mu*(WW1[j][i]-WW2[j][i]+YY[j][i]/mu);
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
			YR[j][t]+=mu*(WW1[j][t]-WW2[j][t]);
		}
	}
	return YR;
}

