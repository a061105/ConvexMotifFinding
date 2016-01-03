#include"incl_head.h"
extern void			ResOut( MAT_D W1);
extern double LossfuncW1(	MAT_D CC,
							MAT_D W1,
							MAT_D W2,
							MAT_D YY);
extern		MAT_D UpdateY(	MAT_D WW1,
							MAT_D WW2,
							MAT_D YY);
extern	MAT_D OptPhaseOne(  MAT_D CC,
							MAT_D WW1,
							MAT_D WW2,
							MAT_D YY,
							int num_k);

int main()
{
	string	    Seq="0110000101110000011100000110110001100101";
	string    InSeq="0110000101110000011100000110110001100101";
	string    word="apple";
	// construct C
	vector<double> col0(Tseq,cost_un),col1(Tseq,0.0);
	MAT_D C(J,col1);
	C[0].assign(col0.begin(),col0.end());
	// initialize W1, all unassigned
	MAT_D W1(C);
	cout<<"Initial is:"<<endl;
	ResOut(W1);
	// continue construct C
	for(int t=0; t<Tseq; t++){
		int dig=0;
		if(InSeq[t]=='1') dig=1;
		for(int kk=0; kk<KG; kk++){
			for(int j=0; j<2*L; j++){
				if(j%2!=dig){	// j odd assign to 0, j even assign to 1
					C[(1+2*L*kk)+j][t]=cost_mis;
				}
			}
		}
	}
	
	// initialize W2, all optimal
	MAT_D W2(J,col1);
	for(int i=1; i<=word_length; i++){ // consider ith char
		int index1=word[i-1]-'a';      // ascii of ith char-97
		for(int j=1; j<=L; j++){	// jth digit of ith char's pattern
			int t=(i-1)*L+j-1;		//pisition in seq
			int dig=0;
			if(Seq[t]=='1') dig=1;	//correct assign should be 1 or 0
			W2[(2*j-2)+dig+(2*L*index1+1)][t]=1;
		}
	}
	cout<<"Optimal is:"<<endl;
	ResOut(W2);
	MAT_D Y(J,col1);  // initialize Y with all zeros
	// start rowlling
	for(int Iter=0; Iter<update_num; Iter++){
		
		for(int Inner_iter=0; Inner_iter<Inner_num; Inner_iter++){
			int k_num=Iter+1;
		W1=OptPhaseOne(C,W1,W2,Y,k_num);
		//cout<<"Inner Step"<<Inner_iter+1<<":"<<endl;
		cout<<"Loss func value:"<<LossfuncW1(C,W1,W2,Y)<<endl;
		}
		Y=UpdateY(W1,W2,Y);
		cout<<"-----End Outer Step"<<Iter+1<<"-----"<<endl;
		cout<<"---Loss func value:"<<LossfuncW1(C,W1,W2,Y)<<"---"<<endl;
		//ResOut(W1);
	}
	cout<<"End output:"<<endl;
	ResOut(W1);
	return 0;
}