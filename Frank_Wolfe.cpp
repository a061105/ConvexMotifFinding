#include"incl_head.h"

vector<int> Frank_Wolfe(	MAT_D Gradf)
{
	vector<int> R(Tseq,MAX_NUMBER);
	vector<int> col_int(Tseq,0);
	vector<double> col_dou(Tseq,MAX_NUMBER);
	vector<vector<int>> B(J,col_int);
	MAT_D M(J,col_dou);
	// initialize M
	M[0][0]=Gradf[0][0];// Un
	for(int i=0; i<KG; i++){   // a1_0,a1_1,b1_0,b1_1...z1_1
		M[2*L*i+1][0]=Gradf[2*L*i+1][0];
		M[2*L*i+2][0]=Gradf[2*L*i+2][0];
	}
	// using DP to find argmin<Gradf,W>
	for(int t=0; t<Tseq-1; t++){		// for all digits
		// for unassigned case
		for(int i=0; i<KG+1; i++){   // for an Un
				int s1=0;
				int s2;
				if(i==0){ 
					s2=0;
					if( M[s2][t+1]>M[s1][t]+Gradf[s2][t+1] ){ // s2 to be Un
						M[s2][t+1]=M[s1][t]+Gradf[s2][t+1];
						B[s2][t+1]=s1;
					}
				}else{
					s2=2*L*(i-1)+1;
					if( M[s2][t+1]>M[s1][t]+Gradf[s2][t+1] ){ // s2 to be a1_0 b1_0...z1_0
						M[s2][t+1]=M[s1][t]+Gradf[s2][t+1];
						B[s2][t+1]=s1;
					}
					s2++;
					if( M[s2][t+1]>M[s1][t]+Gradf[s2][t+1] ){ // s2 to be a1_1 b1_1...z1_1
						M[s2][t+1]=M[s1][t]+Gradf[s2][t+1];
						B[s2][t+1]=s1;
					}
				}
		}
		// for assigned case
		for(int ch=0; ch<KG; ch++){ // for all chars a b c ...
			for(int l=0; l<L-1; l++){ 
					int s1=2*ch*L+1+2*l;
					// for all digit in a char a1_0 a2_0 ... aL-1_0
					for(int s2=s1+2; s2<=s1+3; s2++){				// s2 can only have two choice->  the next one to be 0 or 1
						if( M[s2][t+1]>M[s1][t]+Gradf[s2][t+1] ){  
							M[s2][t+1]=M[s1][t]+Gradf[s2][t+1];
							B[s2][t+1]=s1;
						}
					}
					s1=2*ch*L+2+2*l;
					// for all digit in a char a1_1 a2_1 ... aL-1_1
					for(int s2=s1+1; s2<=s1+2; s2++){				// s2 can only have two choice->  the next one to be 0 or 1
						if( M[s2][t+1]>M[s1][t]+Gradf[s2][t+1] ){  
							M[s2][t+1]=M[s1][t]+Gradf[s2][t+1];
							B[s2][t+1]=s1;
						}
					}

			}
			for(int i=0; i<KG+1; i++){   // for an end char_digit aL(16) bL(32)
				int s1=2*(ch+1)*L;
				int s2;
				if(i==0){
					s2=0;				// s2 can have choice of every U(0)
					if( M[s2][t+1]>M[s1][t]+Gradf[s2][t+1] ){ 
						M[s2][t+1]=M[s1][t]+Gradf[s2][t+1];
						B[s2][t+1]=s1;
					}
				}else{
					s2=2*(i-1)*L+1;		// s2 can have choice of every a1_0(1) a1_1(2) b1(17)...
					if( M[s2][t+1]>M[s1][t]+Gradf[s2][t+1] ){ 
						M[s2][t+1]=M[s1][t]+Gradf[s2][t+1];
						B[s2][t+1]=s1;
					}
					s2++;
					if( M[s2][t+1]>M[s1][t]+Gradf[s2][t+1] ){ 
						M[s2][t+1]=M[s1][t]+Gradf[s2][t+1];
						B[s2][t+1]=s1;
					}
				}
			}
		}
	}
	
	// Generating R: descent direction
	double M_max=MAX_NUMBER;
	for(int i=0; i<J; i++){   // RT=argmin(M(T,S))
		if(M_max>M[i][Tseq-1]){
			R[Tseq-1]=i;
			M_max=M[i][Tseq-1];
		}
	}
	for(int t=Tseq-1; t>=1; t--){
		R[t-1]=B[R[t]][t];
	}
	return R;
}