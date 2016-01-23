#include"incl_head.h"
extern void Report_dir(vector<int> AT);
vector<int> Topk( vector<double> Array, int k);

vector<int> Frank_Wolfe(	MAT_D Gradf)
{
	vector<int> R(Tseq,MAX_NUMBER);
	vector<int> col_int(Tseq,0);
	vector<double> col_dou(Tseq,MAX_NUMBER);
	vector<vector<int> > B(J,col_int);
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
					// for all digit in a char a1_0 a2_0 ... aLmin-1_0
					for(int s2=s1+2; s2<=s1+3; s2++){				// s2 can only have two choice->  the next one to be 0 or 1
						if( M[s2][t+1]>M[s1][t]+Gradf[s2][t+1] ){  
							M[s2][t+1]=M[s1][t]+Gradf[s2][t+1];
							B[s2][t+1]=s1;
						}
					}
					s1=2*ch*L+2+2*l;
					// for all digit in a char a1_1 a2_1 ... aLmin-1_1
					for(int s2=s1+1; s2<=s1+2; s2++){				// s2 can only have two choice->  the next one to be 0 or 1
						if( M[s2][t+1]>M[s1][t]+Gradf[s2][t+1] ){  
							M[s2][t+1]=M[s1][t]+Gradf[s2][t+1];
							B[s2][t+1]=s1;
						}
					}

			}
			for(int l=Lmin-1; l<L; l++){
				for(int i=0; i<KG+1; i++){   // for an possibly end char_digit aLmin...aL bLmin...bL ...
					int s1=2*ch*L+1+2*l;
					int s2;
					if(i==0){
						s2=0;				// s2 can have choice of every U(0)
						if( M[s2][t+1]>M[s1][t]+Gradf[s2][t+1] ){ 
							M[s2][t+1]=M[s1][t]+Gradf[s2][t+1];
							B[s2][t+1]=s1;
						}
					}else{
						s2=2*(i-1)*L+1;		// s2 can have choice of every a1_0(1) a1_1(2) b1(17)...
						if( M[s2][t+1]>=M[s1][t]+Gradf[s2][t+1] ){ 
							M[s2][t+1]=M[s1][t]+Gradf[s2][t+1];
							B[s2][t+1]=s1;
						}
						s2++;
						if( M[s2][t+1]>=M[s1][t]+Gradf[s2][t+1] ){ 
							M[s2][t+1]=M[s1][t]+Gradf[s2][t+1];
							B[s2][t+1]=s1;
						}
					}
					s1++;
					if(i==0){
						s2=0;				// s2 can have choice of every U(0)
						if( M[s2][t+1]>M[s1][t]+Gradf[s2][t+1] ){ 
							M[s2][t+1]=M[s1][t]+Gradf[s2][t+1];
							B[s2][t+1]=s1;
						}
					}else{
						s2=2*(i-1)*L+1;		// s2 can have choice of every a1_0(1) a1_1(2) b1(17)...
						if( M[s2][t+1]>=M[s1][t]+Gradf[s2][t+1] ){ 
							M[s2][t+1]=M[s1][t]+Gradf[s2][t+1];
							B[s2][t+1]=s1;
						}
						s2++;
						if( M[s2][t+1]>=M[s1][t]+Gradf[s2][t+1] ){ 
							M[s2][t+1]=M[s1][t]+Gradf[s2][t+1];
							B[s2][t+1]=s1;
						}
					}
				}
			}
		}
	}
	
	// Generating R: descent direction
	double M_max=M[0][Tseq-1];
	R[Tseq-1]=0;
	for(int kk=0; kk<KG; kk++){
		for(int j=Lmin-1; j<L; j++){   // RT=argmin(M(T,S))
			int i0=2*L*kk+2*j+1;
			if(M_max>M[i0][Tseq-1]){
				R[Tseq-1]=i0;
				M_max=M[i0][Tseq-1];
			}
			int i1=2*L*kk+2*j+2;
			if(M_max>M[i1][Tseq-1]){
				R[Tseq-1]=i1;
				M_max=M[i1][Tseq-1];
			}
		}
	}
	for(int t=Tseq-1; t>=1; t--){
		R[t-1]=B[R[t]][t];
	}
	//Report_dir(R);
	return R;
}


MAT_D GroupConsDir(MAT_D Gradw2, vector<int>& pattern_ind){

	int Lopt=2*Kopt*L+1;
	vector<double> col0(Tseq,0);
	vector<double> pattern_score;
	MAT_D R2(J,col0);
	// Finding a pattern for each pattern slot pa_num=0,1,2 ... Kopt-1
	for(int pa_num=0; pa_num<KG; pa_num++){
		double thispattern=0;
		//cout<<"P"<<pa_num+1<<": ";
		int pa_L=pa_num*2*L+1, pa_R=pa_L+2*L;
		vector<double> pat_length_score;
		for(int dig=0; dig<L; dig++){
			double choose_zero=0, choose_one=0;
			int zero_slot=pa_L+2*dig;
			int one_slot=zero_slot+1;
			vector<double> zero_vec,one_vec;
			for(int t=0; t<Tseq; t++){
				// Calculating score of choosing zero & formulating dir
				if(Gradw2[zero_slot][t]<0){
					choose_zero+=-Gradw2[zero_slot][t];
					zero_vec.push_back(1.0);
				}else{
					zero_vec.push_back(0);
				}
				// Calculating score of choosing one & formulating dir
				if(Gradw2[one_slot][t]<0){
					choose_one+=-Gradw2[one_slot][t];
					one_vec.push_back(1.0);
				}else{
					one_vec.push_back(0);
				}
			}
			if(choose_zero>=choose_one){// then choose zero
				R2[zero_slot].assign(zero_vec.begin(),zero_vec.end());
				thispattern+=choose_zero;
				//cout<<"0";
			}else{						// choose one
				R2[one_slot].assign(one_vec.begin(),one_vec.end());
				thispattern+=choose_one;
				//cout<<"1";
			}
			pat_length_score.push_back(thispattern);
		}// end for this digit
		
		vector<double>::iterator  opt_pos;
		opt_pos=max_element(pat_length_score.begin()+Lmin-1,pat_length_score.end());
		pattern_score.push_back(*opt_pos);
		for(int i=L-1; opt_pos!=(pat_length_score.begin()+i); i--){
			int zero_slot=pa_L+2*i;
			int one_slot=zero_slot+1;
			R2[one_slot]=col0;
			R2[zero_slot]=col0;
			//cout<<'-';
		}
	//cout<<endl;
	}//end for this pattern
	
	pattern_ind=Topk(pattern_score,Kopt);

	return R2;
}


