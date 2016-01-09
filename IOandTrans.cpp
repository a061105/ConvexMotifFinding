#include"incl_head.h"

string Word2Bin(const string word){
	string res;
	for(int p=0; p<word_length; p++){
		int wp=word[p];
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