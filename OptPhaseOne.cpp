#include"incl_head.h"
vector<vector<double>> Gradf(		vector<vector<double>> CC,
					vector<vector<double>> WW1,
					vector<vector<double>> WW2,
					vector<vector<double>> YY);

extern vector<int> Frank_Wolfe(	vector<vector<double>> Gradf);
void Report_dir(vector<int> AT);

vector<vector<double>> OptPhaseOne(	vector<vector<double>> CC,
					vector<vector<double>> WW1,
					vector<vector<double>> WW2,
					vector<vector<double>> YY,
					int num_k)
{
	
	vector<vector<double>> Gradw1;
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

vector<vector<double>> Gradf(	vector<vector<double>> CC,
				vector<vector<double>> WW1,
				vector<vector<double>> WW2,
				vector<vector<double>> YY)
{
	vector<vector<double>> GG=CC;
	for(int i=0; i<Tseq; i++){
		for(int j=0; j<J; j++){
			GG[j][i]+=mu*(WW1[j][i]-WW2[j][i]+YY[j][i]/mu);
		}
	}
	return GG;
}

vector<vector<double>> UpdateY(	vector<vector<double>> WW1,
								vector<vector<double>> WW2,
								vector<vector<double>> YY){
	vector<vector<double>> YR(YY);
	for(int t=0; t<Tseq; t++){
		for(int j=0; j<J; j++){
			YR[j][t]+=mu*(WW1[j][t]-WW2[j][t]);
		}
	}
	return YR;
}

void Report_dir(vector<int> AT){
	cout<<"Current direction is :"<<endl;
	for(int t=0; t<Tseq; t++){
		int Max_stat=AT[t];
		if(Max_stat==0){
			cout<<"Un"<<" "<<endl;
		}else{
			int charnum=(Max_stat-1)/(2*L);
			int inner_stat=Max_stat%(2*L);
			int One_or_zero=(Max_stat%(2*L)+1)%2;  
			char AA='a'+charnum;
			int Inner_seq=(Max_stat%(2*L)+1)/2;
				if(Inner_seq==0) Inner_seq=L;
			char Onezero='0'+One_or_zero;
			cout<<AA<<"_"<<Onezero<<"_"<<Inner_seq<<" "<<endl;
		}
		
	}
	return;
}