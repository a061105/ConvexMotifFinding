#include"incl_head.h"

double OptStep1(MAT_D CC,
				MAT_D W1,
				MAT_D W2,
				MAT_D YY,
				vector<int> Aton){
	double MuBF2=0,CB=0,MuAB=0; // B=S-W1, A=W1-W2+Y/mu, S-descent direction
	
	for(int j=0; j<J; j++){
		for(int t=0; t<Tseq; t++){
			double Sjt=(double)(Aton[t]==j);
			double Bjt=Sjt-W1[j][t];
			double MuAjt=mu*(W1[j][t]-W2[j][t])+YY[j][t];
			MuBF2+=mu*pow(Bjt,2);
			CB+=CC[j][t]*Bjt;
			MuAB+=MuAjt*Bjt;
		}
	}
	if(MuBF2<1e-6) MuBF2=1e-6;
	return -(MuAB+CB)/(MuBF2);
}


double OptStep2(MAT_D W1,
				MAT_D W2,
				MAT_D YY,
				MAT_D DIR,
				vector<int> pattern_ind ){
	double step=0;
	double BF2=0,AB=0; // B=W2-S, A=W1-W2+Y/mu, S-descent direction
		for(int j=0; j<J; j++){
			int kk=(j-1)/(2*L);
			vector<int>::iterator ind=pattern_ind.begin();
			ind=find(pattern_ind.begin(),pattern_ind.end(),kk);
			for(int t=0; t<Tseq; t++){
				double Bjt=W2[j][t];
				if(ind!=pattern_ind.end()) Bjt=-DIR[j][t]+W2[j][t];
				double Ajt=(W1[j][t]-W2[j][t])+YY[j][t]/mu;
				BF2+=pow(Bjt,2);
				AB+=Ajt*Bjt;
			}
		}
		if(BF2<1e-8) BF2=1e-8;
		step=-AB/BF2;
		if(step>1.0) step=1.0;
		if(step<0) step=0;
	return step;
}