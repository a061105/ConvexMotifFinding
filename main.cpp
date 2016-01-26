#include"incl_head.h"
extern void			ResOut( MAT_D W1,
							string InSeq);
extern double LossfuncW1(	MAT_D CC,
							MAT_D W1,
							MAT_D W2,
							MAT_D YY);
extern double LossfuncW1(	MAT_D CC,
							MAT_D W1);
extern		MAT_D UpdateY(	MAT_D WW1,
							MAT_D WW2,
							MAT_D YY);
extern void OptPhaseOne(    MAT_D& CC,
							MAT_D& WW1,
							MAT_D& WW2,
							MAT_D& YY,
							double& step_size);
extern void OptPhaseTwo(	MAT_D& CC,
							MAT_D& WW1,
							MAT_D& WW2,
							MAT_D& YY,
							double& step_size);
extern double diff(			MAT_D WW1,
							MAT_D WW2);
extern double DisToOne(		MAT_D WW1);
extern string Word2Bin(const string word);

int main()
{
	string InSeq=Word2Bin(word);
	// construct C
	vector<double> col0(Tseq,cost_un),col1(Tseq,0.0);
	MAT_D Initial(J,col1);
	MAT_D C(J,col1);
	C[0].assign(col0.begin(),col0.end());
	// initialize W1, all unassigned
	MAT_D W1(C);
	// Initialize W2 -----Unassigned
	MAT_D W2(C);
	
	
	// continue construct C 
	for(int t=0; t<Tseq; t++){
		bool dig=(InSeq[t]=='1');
		for(int kk=0; kk<KG; kk++){
			for(int j=0; j<2*L; j++){
				if(j%2!=dig){	// j odd assign to 0, j even assign to 1
					C[(1+2*L*kk)+j][t]=cost_mis;
				}
			}
		}
	}
	ResOut(W1,InSeq);
	// add prefer for each pattern
	for(int t=0; t<Tseq; t++){
		for(int kk=0; kk<KG; kk++){
			int patnum=kk;
			for(int j=L; j>0; j--){
				int zero_slot=2*L*kk+2*j-1;
				bool ifpref1=patnum%2;
				patnum/=2;
				double rapre=random;
				C[zero_slot+!ifpref1][t]+=prefer;
				if(kk%2==1){
					C[zero_slot][t]+=2*short_prefer;
					C[zero_slot+1][t]+=2*short_prefer;
				}else{
					C[zero_slot][t]+=short_prefer*(0);
					C[zero_slot+1][t]+=short_prefer*(0);
				}
			}
		}
	}
	//----------------------initialize W2, all optimal--------------------------
	/*MAT_D Wopt(J,col1);
	for(int i=1; i<=word_length; i++){ // consider ith char
		int index1=word[i-1]-'a'+1;      // ascii of ith char-97
		index1=index1*pow(2,(double)(L-Lmin));
		for(int j=1; j<=Lmin; j++){	// jth digit of ith char's pattern
			int t=(i-1)*Lmin+j-1;		//pisition in seq
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
		for(int Inner_iter=0; Inner_iter<Inner_num; Inner_iter++){
			double step_length=2.0/(Inner_iter+2);
			OptPhaseOne(C,W1,W2,Y,step_length);
			if(step_length<inner_eps) break;
		}
		//Optimization Phase Two
		W2=Initial;
		for(int Inner_iter=0; Inner_iter<Inner_num; Inner_iter++){
			double step_length=2.0/(Inner_iter+2);
			OptPhaseTwo(C,W1,W2,Y,step_length);
			if(step_length<inner_eps) break;
		}
		//Optimization Phase Three
		Y=UpdateY(W1,W2,Y);
		double diffW12=diff(W1,W2);
		cout<<"L:"<<LossfuncW1(C,W1)<<"_"<<LossfuncW1(C,W2)<<" "<<"D:"<<diffW12<<endl;
		//ResOut(W2);
		if(diffW12<outer_eps) break;
	}
	cout<<"End output W1:"<<endl;
	ResOut(W1,InSeq);
	cout<<"End output W2:"<<endl;
	ResOut(W2,InSeq);
	//cout<<"Optimal loss val:"<<LossfuncW1(C,Wopt)<<endl;
	cout<<"Max entree distance to one:"<<DisToOne(W2)<<endl;
	return 0;
}
