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
extern void SamplePhaseTwo(	MAT_D& CC,
							MAT_D& WW1,
							MAT_D& WW2,
							MAT_D& YY,
							double& step_size,
							vector<int>& Aton,
							double& min_loss);
extern double diff(			MAT_D WW1,
							MAT_D WW2);
extern double DisToOne(		MAT_D WW1);
extern string Word2Bin(const string word);

int main()
{
	string InSeq=Word2Bin(word);
	// construct C
	vector<double> col0(Tseq,cost_un),col1(Tseq,0.0);
	MAT_D Initial(J1+J2,col1);
	MAT_D C(J1+J2,col1);
	C[0].assign(col0.begin(),col0.end());
	// initialize W1, all unassigned
	MAT_D W1(C);
	// Initialize W2 -----Unassigned
	MAT_D W2(C);
	
	
	// continue construct C 
	for(int t=0; t<Tseq; t++){
		bool dig=(InSeq[t]=='1');
		for(int kk=0; kk<KG1; kk++){
			for(int j=0; j<2*L; j++){
				if(j%2!=dig){	// j odd assign to 0, j even assign to 1
					C[(1+2*L*kk)+j][t]=cost_mis;
				}
			}
		}
		for(int kk=0; kk<KG2; kk++){
			for(int j=0; j<2*Lmin; j++){
				if(j%2!=dig){	// j odd assign to 0, j even assign to 1
					C[2*KG1*L+(1+2*Lmin*kk)+j][t]=cost_mis;
				}
			}
		}
	}
	ResOut(W1,InSeq);
	// add prefer for each pattern
	for(int t=0; t<Tseq; t++){
		for(int kk=0; kk<KG1; kk++){
			int patnum=kk;
			C[2*L*kk+1][t]+=2*short_prefer*random;
			C[2*L*kk+2][t]+=2*short_prefer*random;
			for(int j=L; j>0; j--){
				int zero_slot=2*L*kk+2*j-1;
				bool ifpref1=patnum%2;
				patnum/=2;
				C[zero_slot+!ifpref1][t]+=prefer;
			}
		}
		for(int kk=0; kk<KG2; kk++){
			int patnum=kk;
			C[2*Lmin*kk+J1][t]+=2*short_prefer*random;
			C[2*Lmin*kk+J1+1][t]+=2*short_prefer*random;
			for(int j=Lmin; j>0; j--){
				int zero_slot=J1+2*Lmin*kk+2*j-2;//check
				bool ifpref1=patnum%2;
				patnum/=2;
				C[zero_slot+!ifpref1][t]+=prefer;
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
	MAT_D Y(J1+J2,col1);  // initialize Y with all zeros
	double diffW12=1.0;
	// start rowlling
	for(int Iter=0; Iter<update_num; Iter++){
		// Optimization Phase One
		for(int Inner_iter=0; Inner_iter<Inner_num; Inner_iter++){
			double step_length=2.0/(Inner_iter+2);
			OptPhaseOne(C,W1,W2,Y,step_length);
			if(step_length<inner_eps) break;
		}
		//Optimization Phase Two
	//	W2=Initial;
		for(int Inner_iter=0; Inner_iter<Inner_num; Inner_iter++){
			double step_length=2.0/(Inner_iter+2);
			OptPhaseTwo(C,W1,W2,Y,step_length);
			if(step_length<inner_eps) break;
		}
		//Optimization Phase Three
		Y=UpdateY(W1,W2,Y);
		diffW12=diff(W1,W2);
		cout<<"L:"<<LossfuncW1(C,W1)<<"_"<<LossfuncW1(C,W2)<<" "<<"D:"<<diffW12<<endl;
		//ResOut(W2);
		if(diffW12<=outer_eps) break;
	}
	//start sampllig
	if(diffW12>outer_eps){
	W2=Initial;
	vector<int> Aton(Tseq,0);
	double min_loss=MAX_NUMBER;
		for(int Inner_iter=0; Inner_iter<Inner_num*100; Inner_iter++){
			double step_length=1;
			double the_loss=MAX_NUMBER;
			vector<int> InAton(Tseq,0);
			SamplePhaseTwo(C,W1,W2,Y,step_length,InAton,the_loss);
			if(min_loss>the_loss){
				min_loss=the_loss;
				Aton=InAton;
				cout<<"Sample Loss: "<<the_loss<<endl;
			}
			if(min_loss<outer_eps) break;
		}

	MAT_D WS(Initial);
	for(int t=0; t<Tseq; t++){
		int j=Aton[t];
		WS[j][t]+=1;
	}

	/*cout<<"End output W1:"<<endl;
	ResOut(W1,InSeq);
	cout<<"End output W2:"<<endl;
	ResOut(W2,InSeq);*/
		cout<<"End output WS:"<<endl;
		ResOut(WS,InSeq);
	}else{
		cout<<"End output W1"<<endl;
		ResOut(W1,InSeq);
	}

	//cout<<"Optimal loss val:"<<LossfuncW1(C,Wopt)<<endl;
	cout<<"Max entree distance to one:"<<DisToOne(W2)<<endl;
	return 0;
}
