#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <functional>
#include <cstdlib>
#include <list>
#include<algorithm>

using namespace std;
const int Nframe =3;
const int Nrod =200;
const int Nrodp =41;
const int L=140;
const double r=1.1246;
const int imax=100000000; 

int main()
{
    int seqi,seqj,seqk,seqt,seqw,sun;
    int line=0;
    char c;
    char type1,type2;
	
	/// count the lines of the file /////
    fstream fin("outtest_m_n3.dat",ios::in);
    while (fin.get(c))
    {
    	if (c=='\n')
    		line++;
	}
	
    static int cluster[imax][3];
    ifstream infile("outtest_m_n3.dat",ios::in);
	ofstream outfile("cluster.dat",ios::out);
	
	///// read the file and put the data into the array////
	for(seqi=0;seqi<line;seqi++)
	{
		string temp;
		getline(infile,temp);
		stringstream word(temp);
		word>>cluster[seqi][0]>>type1>>cluster[seqi][1]>>type2>>cluster[seqi][2];
	}
	infile.close();
	
	///// deal the array /////
	sun=1; ///count the numbers of cluster //// 
	for(seqj=0;seqj<line;seqj++)
	{	
		list<int> slist;
		if ( cluster[seqj][0]!=cluster[seqj-1][0]){
			sun=1;
		}
		if ( cluster[seqj][1]!=0 ) {
		slist.push_back(cluster[seqj][1]);
		slist.push_back(cluster[seqj][2]);
		list<int>::iterator num;
		for(num=slist.begin();num!=slist.end();num++)  
		{
			seqw=seqj+1;
			while (cluster[seqw][0]==cluster[seqj][0])
			{	
				if (cluster[seqw][1]==*num || cluster[seqw][2]==*num ){
					slist.push_back(cluster[seqw][1]);
					slist.push_back(cluster[seqw][2]);
					cluster[seqw][1]=0;   ///// avoid the list repeat output //////
					cluster[seqw][2]=0;
					seqw++;
				} else {
					seqw++;
				}

			} 
		}
		slist.sort();
		slist.unique();
		list<int>::iterator itor;
		outfile<<cluster[seqj][0]<<","<<slist.size()<<","<<"{ ";
		for(itor=slist.begin();itor!=slist.end();++itor)
			outfile<<*itor<<" ";
		outfile<<"}"<<" All num:"<<sun<<endl;
		slist.clear();
		sun++;
		} else {
			cout<<"see you"<<" \n";
		}
	}
	
	outfile.close();
	return 0;
}
