#include"incl_head.h"

void ResOut( MAT_D W1,
			string InSeq){
	ofstream myfile,logfile;
	vector<string> color;
	string line;
	ifstream colorfile ("color.txt");
	if (colorfile.is_open()){
		while ( getline (colorfile,line) ){
			color.push_back(line);
		}
		colorfile.close();
	}else{
		cout<<"Unable to open file";
	}
	myfile.open ("result.txt");
	logfile.open("word");
	for(int Digit_num=0; Digit_num<Tseq; Digit_num++){
		int Max_stat=0;
		double Max_entree=0;
		for(int Stat=0; Stat<J1+J2;Stat++){
			if(W1[Stat][Digit_num]>Max_entree){
				Max_stat=Stat;
				Max_entree=W1[Stat][Digit_num];
			}
		}
		if(Max_stat==0){
			myfile<<"{\\color{Black}"<<InSeq[Digit_num]<<"}\n";
		}else{
			if(Max_stat<J1){
				int charnum=(Max_stat-1)/(2*L);
				int inner_stat=Max_stat%(2*L);
				int One_or_zero=(Max_stat%(2*L)+1)%2;  
				int AA=1+charnum;
				int Inner_seq=(Max_stat%(2*L)+1)/2;
				if(Inner_seq==0) Inner_seq=L;
				myfile<<"{\\color{"<<color[AA]<<"}"<<One_or_zero<<"}"<<"\n";
				logfile<<"P_"<<AA<<"_"<<Inner_seq<<endl;
			}else{// to be in shorter pattern
				int charnum=(Max_stat-J1)/(2*Lmin)+KG1;
				int inner_stat=(Max_stat-J1+1)%(2*Lmin);
				int One_or_zero=(inner_stat+1)%2;  
				int AA=1+charnum;
				int Inner_seq=(inner_stat+1)/2;
				if(Inner_seq==0) Inner_seq=Lmin;
				myfile<<"{\\color{"<<color[AA]<<"}"<<One_or_zero<<"}"<<"\n";
				logfile<<"P_"<<AA<<"_"<<Inner_seq<<endl;
			}
		}
	}
	cout<<endl;
	myfile.close();
	logfile.close();
	return;
}
	

double LossfuncW1(	MAT_D CC,
					MAT_D W1,
					MAT_D W2,
					MAT_D YY){
	double loss=0;
	for(int j=0; j<J1+J2; j++){
		for(int t=0; t<Tseq; t++){
			loss+=CC[j][t]*W1[j][t]+mu*pow((W1[j][t]-W2[j][t]+YY[j][t]/mu),2);
			//loss-=mu*pow((YY[j][t]/mu),2);
		}
	}
	return loss;
}

// reload lossfunc
double LossfuncW1(	MAT_D CC,
					MAT_D W1){
	double loss=0;
	for(int j=0; j<J1+J2; j++){
		for(int t=0; t<Tseq; t++){
			loss+=CC[j][t]*W1[j][t];;
		}
	}
	return loss;
}

double diff(MAT_D WW1,
			MAT_D WW2){
	double maxdiff=0;
	for(int j=0; j<J1+J2; j++){
		for(int t=0; t<Tseq; t++){
			if(abs(WW1[j][t]-WW2[j][t])>maxdiff){
				maxdiff=abs(WW1[j][t]-WW2[j][t]);
			}
		}
	}
	return maxdiff;
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

double DisToOne(MAT_D WW1){
	double mindiff=1;
	for(int j=0; j<J1+J2; j++){
		for(int t=0; t<Tseq; t++){
			if(abs(WW1[j][t]-1.0)<mindiff){
				mindiff=abs(WW1[j][t]-1.0);
			}
		}
	}
	return mindiff;
}
