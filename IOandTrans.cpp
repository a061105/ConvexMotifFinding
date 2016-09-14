#include"incl_head.h"

string Word2Bin(const string word){
	string res;
	for(int p=0; p<word_length; p++){
		int wp=word[p]-'a'+1;
		for(int l=L-1; l>=0; l--){
			if(wp/(int)pow((double)2,l)){
				res.push_back('1');
				wp-=(int)pow((double)2,l);
			}else{
				res.push_back('0');
			}
		}
	}
	return res;
}

void printW(MAT_D W){
	int k=0;

	while( k<KG){
		cout<<"Slot"<<k<<":"<<endl;
		for(int t=0; t<Tseq; t++){
			for(int j=0; j<L; j++){
				int zero_slot = 2*j+2*k*L+1;
				int one_slot = 2*j+1+2*k*L+1;
				printf("(%.1f,%.1f)",W[zero_slot][t],W[one_slot][t]);
				//cout<<"("<<W[zero_slot][t]<<","<<W[one_slot][t]<<") ";
			}
			cout<<endl;
		}
		k++;
		cout<<endl;
	}
}
