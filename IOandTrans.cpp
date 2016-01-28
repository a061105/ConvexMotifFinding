#include"incl_head.h"

string Word2Bin(const string word){
	string res;
	string DNA;
	int L_char=4;
	for(int p=0; p<word_length; p++){
		int wp=word[p]-'a'+1;
		if(word[p]==' ') wp=(int)pow((double)2,L_char);
		if(word[p]=='v') wp=0;
		for(int l=L_char-1; l>=0; l--){
			if(wp/(int)pow((double)2,l)){
				res.push_back('1');
				DNA.push_back('G');
				wp-=(int)pow((double)2,l);
			}else{
				res.push_back('0');
				DNA.push_back('A');
			}
		}
		if(word[p]==' '){
			res.push_back('1');
			res.push_back('1');
			DNA.push_back('G');
			DNA.push_back('G');
		}
		if(word[p]=='v'){
			res.push_back('0');
			res.push_back('0');
			DNA.push_back('A');
			DNA.push_back('A');
		}
	}
	ofstream str;
	str.open ("InSeqDNA.txt");
	str<<DNA<<endl;
	str.close();
	return res;
}

vector<int> Topk( vector<double> Array, int k){
	vector<int> Top_indi(k,-1);

	for(int indice=0; indice<k; indice++){
		double M_max=-1;
		for(int i=0; i<Array.size(); i++){   // RT=argmin(M(T,S))
			if(M_max<Array[i]){
				Top_indi[indice]=i;
				M_max=Array[i];
			}
		}
		Array[Top_indi[indice]]=-MAX_NUMBER;
	}
	return Top_indi;
}

