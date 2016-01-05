#include"incl_head.h"

void ResOut( MAT_D W1){
	for(int Digit_num=0; Digit_num<Tseq; Digit_num++){
		int Max_stat=0;
		double Max_entree=0;
		for(int Stat=0; Stat<J;Stat++){
			if(W1[Stat][Digit_num]>Max_entree){
				Max_stat=Stat;
				Max_entree=W1[Stat][Digit_num];
			}
		}
		if(Max_stat==0){
			cout<<"Un"<<" "<<endl;
		}else{
			int charnum=(Max_stat-1)/(2*L);
			int inner_stat=Max_stat%(2*L);
			int One_or_zero=(Max_stat%(2*L)+1)%2;  
			int AA=1+charnum;
			int Inner_seq=(Max_stat%(2*L)+1)/2;
				if(Inner_seq==0) Inner_seq=L;
			cout<<"P"<<AA<<"_"<<One_or_zero<<"_"<<Inner_seq<<endl;
		}
	}
	cout<<endl;
	return;
}
	

double LossfuncW1(	MAT_D CC,
					MAT_D W1,
					MAT_D W2,
					MAT_D YY){
	double loss=0;
	for(int t=0; t<Tseq; t++){
		for(int j=0; j<J; j++){
			loss+=CC[j][t]*W1[j][t]+mu*pow((W1[j][t]-W2[j][t]+YY[j][t]/mu),2);
		}
	}
	return loss;
}
					

// Report descent direction after Franke-Wolfe
void Report_dir(vector<int> AT){
	cout<<"Current direction is :"<<endl;
	for(int t=0; t<Tseq; t++){
		int Max_stat=AT[t];
		if(Max_stat==0){
			cout<<"Un"<<endl;
		}else{
			int charnum=(Max_stat-1)/(2*L);
			int inner_stat=Max_stat%(2*L);
			int One_or_zero=(Max_stat%(2*L)+1)%2;  
			int AA=1+charnum;
			int Inner_seq=(Max_stat%(2*L)+1)/2;
				if(Inner_seq==0) Inner_seq=L;
			cout<<"p"<<AA<<"_"<<One_or_zero<<"_"<<Inner_seq<<endl;
		}
		
	}
	return;
}			

