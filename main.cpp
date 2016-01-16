#include"incl_head.h"
extern void			ResOut( MAT_D W1);
extern double LossfuncW1(	MAT_D CC,
							MAT_D W1,
							MAT_D W2,
							MAT_D YY);
extern double LossfuncW1(	MAT_D CC,
							MAT_D W1);
extern		MAT_D UpdateY(	MAT_D WW1,
							MAT_D WW2,
							MAT_D YY);
extern	MAT_D OptPhaseOne(  MAT_D CC,
							MAT_D WW1,
							MAT_D WW2,
							MAT_D YY,
							int num_k);
extern MAT_D OptPhaseTwo(	MAT_D CC,
							MAT_D WW1,
							MAT_D WW2,
							MAT_D YY,
							int num_k);
extern double diff(			MAT_D WW1,
							MAT_D WW2);
extern double DisToOne(		MAT_D WW1);
extern string Word2Bin(const string word);

int main()
{
	string InSeq="011001111010111010011001";
	// construct C
	vector<double> col0(Tseq,cost_un),col1(Tseq,0.0);
	MAT_D C(J,col1);
	C[0].assign(col0.begin(),col0.end());
	// initialize W1, all unassigned
	MAT_D W1(C);
	cout<<"Initial is:"<<endl;
	ResOut(W1);
	// Initialize W2 -----Unassigned
	MAT_D W2(C);
	//adding a little random bias on patterns
	vector<double> Errcol(KG,0.0);
	MAT_D adderr(word_length,Errcol);
	double eps=0.00/KG/word_length;
	srand((unsigned)time(0));
	for(int k=0; k<KG; k++){
		for(int cha=0; cha<word_length; cha++){
			//adderr[cha][k]=eps*random;
			adderr[cha][k]+=2*eps*cha;
		}
	}
	
	
	// continue construct C 
	for(int t=0; t<Tseq; t++){
		int dig=0;
		int cha=t/L;
		if(InSeq[t]=='1') dig=1;
		for(int kk=0; kk<KG; kk++){
			for(int j=0; j<2*L; j++){
				if(j%2!=dig){	// j odd assign to 0, j even assign to 1
					C[(1+2*L*kk)+j][t]=cost_mis+adderr[cha][kk];
				}else{
					C[(1+2*L*kk)+j][t]=adderr[cha][kk];
				}
			}
		}
	}
	// add prefer for each pattern
	for(int t=0; t<Tseq; t++){
		for(int kk=0; kk<KG; kk++){
			int patnum=kk+1;
			for(int j=L; j>0; j--){
				int zero_slot=2*L*kk+2*j-1;
				bool ifpref1=patnum%2;
				patnum/=2;
				C[zero_slot+!ifpref1][t]+=prefer;
			}
		}
	}

	//----------------------initialize W2, all optimal--------------------------
	/*MAT_D Wopt(J,col1);
	for(int i=1; i<=word_length; i++){ // consider ith char
		int index1=word[i-1]-'a';      // ascii of ith char-97
		for(int j=1; j<=L; j++){	// jth digit of ith char's pattern
			int t=(i-1)*L+j-1;		//pisition in seq
			int dig=0;
			if(InSeq[t]=='1') dig=1;	//correct assign should be 1 or 0
			Wopt[(2*j-2)+dig+(2*L*index1+1)][t]=1;
		}
	}
	double Opt_val=LossfuncW1(C,Wopt);
	cout<<"Optimal is:"<<Opt_val<<endl;*/
	//----------------------------------------------------------------------------
	MAT_D Y(J,col1);  // initialize Y with all zeros
	// start rowlling
	for(int Iter=0; Iter<update_num; Iter++){
		// Optimization Phase One
		//cout<<"Phase One"<<endl;
		for(int Inner_iter=0; Inner_iter<Inner_num; Inner_iter++){
			int k_num=Iter+1;
		W1=OptPhaseOne(C,W1,W2,Y,k_num);
		}
		//Optimization Phase Two
		//cout<<"Phase two"<<endl;
		for(int Inner_iter=0; Inner_iter<Inner_num; Inner_iter++){
			int k_num=Iter+1;
		W2=OptPhaseTwo(C,W1,W2,Y,Inner_iter);
		}
		//Optimization Phase Three
		Y=UpdateY(W1,W2,Y);
		//cout<<"-----End Outer Step"<<Iter+1<<"-----"<<endl;
		double diffW12=diff(W1,W2);
		cout<<"L:"<<LossfuncW1(C,W1)<<"_"<<LossfuncW1(C,W2)<<" "<<"D:"<<diffW12<<endl;
		//ResOut(W2);
		if(diffW12<1e-5) break;
	}
	cout<<"End output W1:"<<endl;
	ResOut(W1);
	cout<<"End output W2:"<<endl;
	ResOut(W2);
	//cout<<"Optimal loss val:"<<LossfuncW1(C,Wopt)<<endl;
	cout<<"Max entree distance to one:"<<DisToOne(W1)<<endl;
	return 0;
}