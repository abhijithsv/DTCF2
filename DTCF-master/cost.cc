#include "cost.h"


using namespace std;


cost::cost(Tensor* &Ainit, Tensor* &Binit, Tensor* &Cinit, Grid* &Ginit)
    {
        A=Ainit;
        B=Binit;
        C=Cinit;
	G=Ginit;
        dimA=A->dims;
        dimB=B->dims;
        dimC=C->dims;
        dimG=G->grid_dims;
	nprocs=G->nprocs;
	found=false;
	
    }


void cost::bestCost(int* &idmapA, int* &idmapB, int* &idmapC, int* &ndimG, int* &npgrid)
{	
	
	ndimG = new int[1];
	if(dimA==2 && dimB==2 && dimC==2)
		contract_2D_2D_2D(idmapA, idmapB, idmapC, ndimG, npgrid);
	else if(dimA==4 && dimB==4 && dimC==4)
		contract_4D_4D_4D(idmapA, idmapB, idmapC, ndimG, npgrid);
	else if(dimA==4 && dimB==2 && dimC==4)
		contract_4D_2D_4D(idmapA, idmapB, idmapC, ndimG, npgrid);
	else if(dimA==4 && dimB==4 && dimC==2)
		contract_4D_4D_2D(idmapA, idmapB, idmapC, ndimG, npgrid);		
	else if(dimA==2 && dimB==4 && dimC==4)
		contract_2D_4D_4D(idmapA, idmapB, idmapC, ndimG, npgrid);

 

}

void cost::contract_2D_2D_2D(int* &idmapA, int* &idmapB, int* &idmapC, int* &ndimG, int* &npgrid)
{	
	
	idmapA = new int[2];
	idmapB=new int[2];
	idmapC=new int[2];
	ndimG[0]=2; //Always 2D grid for these contractions
	npgrid=new int[ndimG[0]];
	//Assuming nprocs is a perfect square
	for(int i=0; i<ndimG[0]; i++)	npgrid[i]=sqrt(nprocs);

	//Find where to distribute the contraction indices
	for(int i=0; i<dimA; i++)
	    	{
        		for(int j=0; j<dimB; j++)
       		 	{
        	    		if(A->contr_dim_str[i].compare(B->contr_dim_str[j]) == 0)
          			{
               		 		// i in A and j in B are contracting indices
                			idmapA[i]=0;	idmapA[j]=1;	
					idmapB[j]=1;	idmapB[j]=0;
					idmapC[i]=0;	idmapC[j]=1;	
            			}
        		}
    		}

		
}

void cost::contract_4D_4D_4D(int* &idmapA, int* &idmapB, int* &idmapC, int* &ndimG, int* &npgrid)
{	
	
	ndimG[0]=4; //Always 4D grid for these contractions
	npgrid=new int[ndimG[0]];
	
	int fsize=2*ceil(sqrt(nprocs));
	int* temp = new int[fsize];
	int nfacts=0;
	

	for(int i=1; i<ceil(sqrt(nprocs));i++) if(nprocs%i==0) {temp[nfacts]=i; nfacts++; temp[nfacts]=nprocs/i; nfacts++;}
	
	int *facts = new int[4*nfacts];
	for(int i=0;i<nfacts;i++){for(int j=0; j<4; j++)facts[i*4+j]=temp[i];}
	
	int* curgrid = new int[4];
	found=false;
	
	
	Grid_Comb_4D(facts,curgrid,0,4*nfacts-1,0,4, idmapA, idmapB, idmapC, npgrid);
	
}


void cost::Grid_Comb_4D(int arr[], int data[], int start, int end, int index, int r,int* &idmapA, int* &idmapB, int* &idmapC, int* &npgrid)
{	
	if(found)return;
	if(index==r)
	{	
		int temp=1;
		for (int i=0;i<4;i++) temp=temp*data[i];
		if(temp==nprocs && data[0]==data[1] && data[2]==data[3]) //Additional conditions to ensure symmetric dimensions are in same sized dimensions
		{	
			best_Cost_GGrid_4D(data, idmapA, idmapB, idmapC);
				     	
			for(int i=0;i<4;i++) npgrid[i]=data[i];
			found=true;
		}
	return;
	}

	for (int i=start; i<=end && end-i+1 >= r-index; i++)
    	{
        	data[index] = arr[i];
        	Grid_Comb_4D(arr, data, i+1, end, index+1, r, idmapA,idmapB,idmapC,npgrid);
  
        	// Remove duplicates
        	while (arr[i] == arr[i+1])
             	i++;
    	}

}


void cost::best_Cost_GGrid_4D(int curGrid[], int* &curidmapA, int* &curidmapB, int* &curidmapC)
{	
	curidmapA = new int[4];
	curidmapB=new int[4];
	curidmapC=new int[4];
	int contA[2]={};
	int extA[2]={};
	int contB[2]={};
	int extB[2]={};

int ccnt=0, ecnt=0;
	for(int i=0; i<dimA; i++)
	    	{
        		for(int j=0; j<dimB; j++)
       		 	{
        	    		if(A->contr_dim_str[i].compare(B->contr_dim_str[j]) == 0)
          			{
               		 		// i in A and j in B are contracting indices
                			contA[ccnt]=i; contB[ccnt]=j;
					ccnt++;
						
            			}
        		}
    		}

for(int i=0;i<4;i++) {if(i!=contA[0] && i!=contA[1]) extA[ecnt]=i; ecnt++;}
ecnt=0;
for(int i=0;i<4;i++) {if(i!=contB[0] && i!=contB[1]) extB[ecnt]=i; ecnt++;}

int lar=1;
if(curGrid[3] > curGrid[1]) lar=3;

int size_A=1, size_B=1;
for(int i=0;i<4;i++){ size_A*=A->tensor_size[i];	size_B*=B->tensor_size[i];}

if(size_A > size_B)
{	cout<<"1"<<endl;
	
	if(A->tensor_size[extA[0]] >= A->tensor_size[extA[1]])
	{cout<<"2"<<endl;
		curidmapA[extA[0]]=lar;
		if(A->tensor_str.at(extA[0]) == 'c' || A->tensor_str.at(extA[0]) == 'd' || A->tensor_str.at(extA[0]) == A->tensor_str.at(extA[1]))
			{ cout<<"3"<<endl;
			if(A->tensor_str.at(extA[1]) == 'c' || A->tensor_str.at(extA[1]) == 'd')
				{cout<<"4"<<endl;
				curidmapA[extA[1]]=lar-1; curidmapA[contA[0]]=4-lar; curidmapA[contA[1]] = 4-lar-1;
			  	curidmapB[extB[0]]=4-lar; curidmapB[extB[1]]=4-lar-1; curidmapB[contB[0]]=lar; curidmapB[contB[1]]=lar-1;
				curidmapC[0]=lar; curidmapC[1]=lar-1; curidmapC[2]=4-lar; curidmapC[3]=4-lar-1; 
				}
			else if((A->tensor_str.at(extA[1]) == 'a' && A->tensor_str.at(contA[0]) == 'A') || (A->tensor_str.at(extA[1]) == 'b' && A->tensor_str.at(contA[0]) == 'B'))
				{cout<<"1"<<endl;
				curidmapA[extA[1]]=4-lar; curidmapA[contA[0]]=4-lar-1; curidmapA[contA[1]] =lar-1;
				  curidmapB[extB[0]]=4-lar-1; curidmapB[extB[1]]=lar-1; curidmapB[contB[0]]=4-lar; curidmapB[contB[1]]=lar;
				curidmapC[0]=lar; curidmapC[1]=4-lar; curidmapC[2]=4-lar-1; curidmapC[3]=lar-1; 
				}
			else if((A->tensor_str.at(extA[1]) == 'a' && A->tensor_str.at(contA[1]) == 'A') || (A->tensor_str.at(extA[1]) == 'b' && A->tensor_str.at(contA[1]) == 'B'))
				{cout<<"1"<<endl;
				curidmapA[extA[1]]=4-lar; curidmapA[contA[1]]=4-lar-1; curidmapA[contA[0]] =lar-1;
				  curidmapB[extB[0]]=4-lar-1; curidmapB[extB[1]]=lar-1; curidmapB[contB[0]]=4-lar; curidmapB[contB[1]]=lar;
				curidmapC[0]=lar; curidmapC[1]=4-lar; curidmapC[2]=4-lar-1; curidmapC[3]=lar-1; 
				}
			}
		else if((A->tensor_str.at(extA[0]) == 'a' && A->tensor_str.at(contA[0]) == 'A') || (A->tensor_str.at(extA[0]) == 'b' && A->tensor_str.at(contA[0]) == 'B'))
			{ cout<<"1"<<endl;
				curidmapA[contA[0]]=lar-1; curidmapA[contA[1]]=4-lar; curidmapA[extA[1]]=4-lar-1;
			  curidmapB[extB[0]]=4-lar; curidmapB[extB[1]]=lar-1; curidmapB[contB[0]]=4-lar-1; curidmapB[contB[1]]=lar;
				curidmapC[0]=lar; curidmapC[1]=4-lar-1; curidmapC[2]=4-lar; curidmapC[3]=lar-1; 
			}
		else if((A->tensor_str.at(extA[0]) == 'a' && A->tensor_str.at(contA[1]) == 'A') || (A->tensor_str.at(extA[0]) == 'b' && A->tensor_str.at(contA[1]) == 'B')) 
			{ cout<<"1"<<endl;
			curidmapA[contA[1]]=lar-1; curidmapA[contA[0]]=4-lar; curidmapA[extA[1]]=4-lar-1;
			  curidmapB[extB[0]]=4-lar; curidmapB[extB[1]]=lar-1; curidmapB[contB[0]]=4-lar-1; curidmapB[contB[1]]=lar;
				curidmapC[0]=lar; curidmapC[1]=4-lar-1; curidmapC[2]=4-lar; curidmapC[3]=lar-1; 
			}
	}
	else
	{cout<<"11"<<endl;
		curidmapA[extA[1]]=lar;
		if(A->tensor_str.at(extA[1]) == 'c' || A->tensor_str.at(extA[1]) == 'd' || A->tensor_str.at(extA[0]) == A->tensor_str.at(extA[1]))
			{ 
			if(A->tensor_str.at(extA[0]) == 'c' || A->tensor_str.at(extA[0]) == 'd')
				{curidmapA[extA[0]]=lar-1; curidmapA[contA[0]]=4-lar; curidmapA[contA[1]] = 4-lar-1;
			  	curidmapB[extB[0]]=4-lar; curidmapB[extB[1]]=4-lar-1; curidmapB[contB[0]]=lar; curidmapB[contB[1]]=lar-1;
				curidmapC[0]=lar-1; curidmapC[1]=lar; curidmapC[2]=4-lar; curidmapC[3]=4-lar-1; 
				}
			else if((A->tensor_str.at(extA[0]) == 'a' && A->tensor_str.at(contA[0]) == 'A') || (A->tensor_str.at(extA[0]) == 'b' && A->tensor_str.at(contA[0]) == 'B'))
				{curidmapA[extA[0]]=4-lar; curidmapA[contA[0]]=4-lar-1; curidmapA[contA[1]] =lar-1;
				  curidmapB[extB[0]]=4-lar-1; curidmapB[extB[1]]=lar-1; curidmapB[contB[0]]=4-lar; curidmapB[contB[1]]=lar;
				curidmapC[0]=4-lar; curidmapC[1]=lar; curidmapC[2]=4-lar-1; curidmapC[3]=lar-1; 
				}
			else if((A->tensor_str.at(extA[0]) == 'a' && A->tensor_str.at(contA[1]) == 'A') || (A->tensor_str.at(extA[0]) == 'b' && A->tensor_str.at(contA[1]) == 'B'))
				{curidmapA[extA[0]]=4-lar; curidmapA[contA[1]]=4-lar-1; curidmapA[contA[0]] =lar-1;
				  curidmapB[extB[0]]=4-lar-1; curidmapB[extB[1]]=lar-1; curidmapB[contB[0]]=4-lar; curidmapB[contB[1]]=lar;
				curidmapC[0]=4-lar; curidmapC[1]=lar; curidmapC[2]=4-lar-1; curidmapC[3]=lar-1; 
				}
			}
		else if((A->tensor_str.at(extA[1]) == 'a' && A->tensor_str.at(contA[0]) == 'A') || (A->tensor_str.at(extA[1]) == 'b' && A->tensor_str.at(contA[0]) == 'B'))
			{ curidmapA[contA[0]]=lar-1; curidmapA[contA[1]]=4-lar; curidmapA[extA[0]]=4-lar-1;
			  curidmapB[extB[0]]=4-lar; curidmapB[extB[1]]=lar-1; curidmapB[contB[0]]=4-lar-1; curidmapB[contB[1]]=lar;
				curidmapC[0]=4-lar-1; curidmapC[1]=lar; curidmapC[2]=4-lar; curidmapC[3]=lar-1; 
			}
		else if((A->tensor_str.at(extA[1]) == 'a' && A->tensor_str.at(contA[1]) == 'A') || (A->tensor_str.at(extA[1]) == 'b' && A->tensor_str.at(contA[1]) == 'B')) 
			{ curidmapA[contA[1]]=lar-1; curidmapA[contA[0]]=4-lar; curidmapA[extA[0]]=4-lar-1;
			  curidmapB[extB[0]]=4-lar; curidmapB[extB[1]]=lar-1; curidmapB[contB[0]]=4-lar-1; curidmapB[contB[1]]=lar;
				curidmapC[0]=4-lar-1; curidmapC[1]=lar; curidmapC[2]=4-lar; curidmapC[3]=lar-1; 
			}
	}	
}	
else
{
	if(B->tensor_size[extB[0]] >= B->tensor_size[extB[1]])
	{cout<<"133"<<endl;
		curidmapB[extB[0]]=lar;
		if(B->tensor_str.at(extB[0]) == 'c' || B->tensor_str.at(extB[0]) == 'd' || B->tensor_str.at(extB[0]) == B->tensor_str.at(extB[1]))
			{ cout<<"1313"<<endl;
			if(B->tensor_str.at(extB[1]) == 'c' || B->tensor_str.at(extB[1]) == 'd')
				{curidmapB[extB[1]]=lar-1; curidmapB[contB[0]]=4-lar; curidmapB[contB[1]] = 4-lar-1;
			  	curidmapA[extA[0]]=4-lar; curidmapA[extA[1]]=4-lar-1; curidmapA[contA[0]]=lar; curidmapA[contA[1]]=lar-1;
				curidmapC[0]=4-lar; curidmapC[1]=4-lar-1; curidmapC[2]=lar; curidmapC[3]=lar-1; 
			}
			else if((B->tensor_str.at(extB[1]) == 'a' && B->tensor_str.at(contB[0]) == 'A') || (B->tensor_str.at(extB[1]) == 'b' && B->tensor_str.at(contB[0]) == 'B'))
				{cout<<"1233"<<endl;
				curidmapB[extB[1]]=4-lar; curidmapB[contB[0]]=4-lar-1; curidmapB[contB[1]] =lar-1;
				  curidmapA[extA[0]]=4-lar-1; curidmapA[extA[1]]=lar-1; curidmapA[contA[0]]=4-lar; curidmapA[contA[1]]=lar;
				curidmapC[0]=4-lar-1; curidmapC[1]=lar-1; curidmapC[2]=lar; curidmapC[3]=4-lar; 
				}
			else if((B->tensor_str.at(extB[1]) == 'a' && B->tensor_str.at(contB[1]) == 'A') || (B->tensor_str.at(extB[1]) == 'b' && B->tensor_str.at(contB[1]) == 'B'))
				{cout<<"1333"<<endl;
				curidmapB[extB[1]]=4-lar; curidmapB[contB[1]]=4-lar-1; curidmapB[contB[0]] =lar-1;
				  curidmapA[extA[0]]=4-lar-1; curidmapA[extA[1]]=lar-1; curidmapA[contA[0]]=4-lar; curidmapA[contA[1]]=lar;
				curidmapC[0]=4-lar-1; curidmapC[1]=lar-1; curidmapC[2]=lar; curidmapC[3]=4-lar; 
				}
			}
		else if((B->tensor_str.at(extB[0]) == 'a' && B->tensor_str.at(contB[0]) == 'A') || (B->tensor_str.at(extB[0]) == 'b' && B->tensor_str.at(contB[0]) == 'B'))
			{ cout<<"1343"<<endl;
			curidmapB[contB[0]]=lar-1; curidmapB[contB[1]]=4-lar; curidmapB[extB[1]]=4-lar-1;
			  curidmapA[extA[0]]=4-lar; curidmapA[extA[1]]=lar-1; curidmapA[contA[0]]=4-lar-1; curidmapA[contA[1]]=lar;
				curidmapC[0]=4-lar; curidmapC[1]=lar-1; curidmapC[2]=lar; curidmapC[3]=4-lar-1; 
			}
		else if((B->tensor_str.at(extB[0]) == 'a' && B->tensor_str.at(contB[1]) == 'A') || (B->tensor_str.at(extB[0]) == 'b' && B->tensor_str.at(contB[1]) == 'B'))
			{ cout<<"1353"<<endl;
			curidmapB[contB[1]]=lar-1; curidmapB[contB[0]]=4-lar; curidmapB[extB[1]]=4-lar-1;
			  curidmapA[extA[0]]=4-lar; curidmapA[extA[1]]=lar-1; curidmapA[contA[0]]=4-lar-1; curidmapA[contA[1]]=lar;
				curidmapC[0]=4-lar; curidmapC[1]=lar-1; curidmapC[2]=lar; curidmapC[3]=4-lar-1; 
			}
	}
	else
	{cout<<"155"<<endl;
		curidmapB[extB[1]]=lar;
		if(B->tensor_str.at(extB[1]) == 'c' || B->tensor_str.at(extB[1]) == 'd' || B->tensor_str.at(extB[0]) == B->tensor_str.at(extB[1]))
			{ 
			if(B->tensor_str.at(extB[0]) == 'c' || B->tensor_str.at(extB[0]) == 'd')
				{curidmapB[extB[0]]=lar-1; curidmapB[contB[0]]=4-lar; curidmapB[contB[1]] = 4-lar-1;
			  	curidmapA[extA[0]]=4-lar; curidmapA[extA[1]]=4-lar-1; curidmapA[contA[0]]=lar; curidmapA[contA[1]]=lar-1;
				curidmapC[0]=4-lar; curidmapC[1]=4-lar-1; curidmapC[2]=lar-1; curidmapC[3]=lar; 
				}
			else if((B->tensor_str.at(extB[0]) == 'a' && B->tensor_str.at(contB[0]) == 'A') || (B->tensor_str.at(extB[0]) == 'b' && B->tensor_str.at(contB[0]) == 'B'))
				{curidmapB[extB[0]]=4-lar; curidmapB[contB[0]]=4-lar-1; curidmapB[contB[1]] =lar-1;
				  curidmapA[extA[0]]=4-lar-1; curidmapA[extA[1]]=lar-1; curidmapA[contA[0]]=4-lar; curidmapA[contA[1]]=lar;
				curidmapC[0]=4-lar-1; curidmapC[1]=lar-1; curidmapC[2]=4-lar; curidmapC[3]=lar; 
				}
			else if((B->tensor_str.at(extB[0]) == 'a' && B->tensor_str.at(contB[1]) == 'A') || (B->tensor_str.at(extB[0]) == 'b' && B->tensor_str.at(contB[1]) == 'B'))
				{curidmapB[extB[0]]=4-lar; curidmapB[contB[1]]=4-lar-1; curidmapB[contB[0]] =lar-1;
				  curidmapA[extA[0]]=4-lar-1; curidmapA[extA[1]]=lar-1; curidmapA[contA[0]]=4-lar; curidmapA[contA[1]]=lar;
				curidmapC[0]=4-lar-1; curidmapC[1]=lar-1; curidmapC[2]=4-lar; curidmapC[3]=lar; 
				}
			}
		else if((B->tensor_str.at(extB[1]) == 'a' && B->tensor_str.at(contB[0]) == 'A') || (B->tensor_str.at(extB[1]) == 'b' && B->tensor_str.at(contB[0]) == 'B'))
			{ curidmapB[contB[0]]=lar-1; curidmapB[contB[1]]=4-lar; curidmapB[extB[0]]=4-lar-1;
			  curidmapA[extA[0]]=4-lar; curidmapA[extA[1]]=lar-1; curidmapA[contA[0]]=4-lar-1; curidmapA[contA[1]]=lar;
				curidmapC[0]=4-lar; curidmapC[1]=lar-1; curidmapC[2]=4-lar-1; curidmapC[3]=lar; 
			}
		else if((B->tensor_str.at(extB[1]) == 'a' && B->tensor_str.at(contB[1]) == 'A') || (B->tensor_str.at(extB[1]) == 'b' && B->tensor_str.at(contB[1]) == 'B')) 
			{ curidmapB[contB[1]]=lar-1; curidmapB[contB[0]]=4-lar; curidmapB[extB[0]]=4-lar-1;
			  curidmapA[extA[0]]=4-lar; curidmapA[extA[1]]=lar-1; curidmapA[contA[0]]=4-lar-1; curidmapA[contA[1]]=lar;
				curidmapC[0]=4-lar; curidmapC[1]=lar-1; curidmapC[2]=4-lar-1; curidmapC[3]=lar; 
			}
	}

}	
}




void cost::contract_4D_2D_4D(int* &idmapA, int* &idmapB, int* &idmapC, int* &ndimG, int* &npgrid)
{	
	
	ndimG[0]=4; //Always 4D grid for these contractions
	npgrid=new int[ndimG[0]];
	
	int fsize=2*ceil(sqrt(nprocs));
	int* temp = new int[fsize];
	int nfacts=0;

	for(int i=1; i<ceil(sqrt(nprocs));i++) if(nprocs%i==0) {temp[nfacts]=i; nfacts++; temp[nfacts]=nprocs/i; nfacts++;}
	int *facts = new int[4*nfacts];
	for(int i=0;i<nfacts;i++){for(int j=0; j<4; j++)facts[i*4+j]=temp[i];}
	
	int* curgrid = new int[4];
	found=false;
	Grid_Comb_4Db(facts,curgrid,0,4*nfacts-1,0,4, idmapA, idmapB, idmapC, npgrid);
}



void cost::Grid_Comb_4Db(int arr[], int data[], int start, int end, int index, int r,int* &idmapA, int* &idmapB, int* &idmapC, int* &npgrid)
{	
	if(found)return;
	if(index==r)
	{
		int temp=1;
		for (int i=0;i<4;i++) temp=temp*data[i];
		if(temp==nprocs && data[0]==data[1] && data[2]==data[3]) //Additional conditions to ensure symmetric dimensions are in same sized dimensions
		{	
			best_Cost_GGrid_4Db(data, idmapA, idmapB, idmapC);	     	
			for(int i=0;i<4;i++) npgrid[i]=data[i];
			found=true;
		}
	return;
	}

	for (int i=start; i<=end && end-i+1 >= r-index; i++)
    	{
        	data[index] = arr[i];
        	Grid_Comb_4Db(arr, data, i+1, end, index+1, r, idmapA, idmapB, idmapC, npgrid);
  
        	// Remove duplicates
        	while (arr[i] == arr[i+1])
             	i++;
    	}

}


void cost::best_Cost_GGrid_4Db(int* &curGrid, int* &curidmapA, int* &curidmapB, int* &curidmapC)
{	
	
	curidmapA = new int[4];
	curidmapB=new int[2];
	curidmapC=new int[4];
	int contA[1]={};
	int extA[3]={};
	int contB[1]={};
	int extB[1]={};

	int ccnt=0, ecnt=0;
	for(int i=0; i<dimA; i++)
	    	{
        		for(int j=0; j<dimB; j++)
       		 	{
        	    		if(A->contr_dim_str[i].compare(B->contr_dim_str[j]) == 0)
          			{
               		 		// i in A and j in B are contracting indices
                			contA[ccnt]=i; contB[ccnt]=j;
					ccnt++;
						
            			}
        		}
    		}

	for(int i=0;i<4;i++) {if(i!=contA[0]) extA[ecnt]=i; ecnt++;}
	ecnt=0;
	for(int i=0;i<4;i++) {if(i!=contB[0]) extB[ecnt]=i; ecnt++;}

	int lar=1;
	if(curGrid[3] > curGrid[1]) lar=3;

	int size_A=1; int size_B=1;
	for(int i=0;i<4;i++){ size_A*=A->tensor_size[i];	size_B*=B->tensor_size[i];}
	
	int temp;
	for(int i=1; i<3; i++){if(A->tensor_size[extA[i]] > A->tensor_size[extA[0]] ) {temp = extA[0]; extA[0]=extA[i]; extA[i]=temp;}}
	
	
	if(A->tensor_size[extA[0]] >= A->tensor_size[extA[1]] && A->tensor_size[extA[0]] >= A->tensor_size[extA[2]])
	{
		curidmapA[extA[0]]=lar;
		if(A->tensor_str.at(extA[0]) == 'c' || A->tensor_str.at(extA[0]) == 'd')
		{ 
			if(A->tensor_size[extA[1]] >= A->tensor_size[extA[2]])
			{
				if(A->tensor_str.at(extA[1]) == 'c' || A->tensor_str.at(extA[1]) == 'd')
				{
					curidmapA[extA[1]]=lar-1; curidmapA[extA[2]]=4-lar; curidmapA[contA[0]]=4-lar-1;
					curidmapB[contB[0]]=4-lar; curidmapB[extB[0]]=4-lar-1;
					curidmapC[0]=lar; curidmapC[1]=lar-1; curidmapC[2]=4-lar; curidmapC[3]=4-lar-1;
				}
				else if(A->tensor_str.at(extA[1]) == A->tensor_str.at(extA[2]))
				{
					curidmapA[extA[1]]=4-lar; curidmapA[extA[2]]=4-lar-1; curidmapA[contA[0]]=lar-1;
					curidmapB[contB[0]]=lar; curidmapB[extB[0]]=lar-1;
					curidmapC[0]=lar; curidmapC[1]=4-lar; curidmapC[2]=4-lar-1; curidmapC[3]=lar-1;	 
				}
				else
				{
					curidmapA[extA[1]]=4-lar; curidmapA[extA[2]]=lar-1; curidmapA[contA[0]]=4-lar-1;
					curidmapB[contB[0]]=4-lar; curidmapB[extB[0]]=4-lar-1;
					curidmapC[0]=lar; curidmapC[1]=4-lar; curidmapC[2]=lar-1; curidmapC[3]=4-lar-1;
				}
			}
			else if(A->tensor_size[extA[2]] >= A->tensor_size[extA[1]])
			{
				if(A->tensor_str.at(extA[2]) == 'c' || A->tensor_str.at(extA[2]) == 'd')
				{
					curidmapA[extA[2]]=lar-1; curidmapA[extA[1]]=4-lar; curidmapA[contA[0]]=4-lar-1;
					curidmapB[contB[0]]=4-lar; curidmapB[extB[0]]=4-lar-1;
					curidmapC[0]=lar; curidmapC[2]=lar-1; curidmapC[1]=4-lar; curidmapC[3]=4-lar-1;
				}
				else if(A->tensor_str.at(extA[2]) == A->tensor_str.at(extA[1]))
				{
					curidmapA[extA[2]]=4-lar; curidmapA[extA[1]]=4-lar-1; curidmapA[contA[0]]=lar-1;
					curidmapB[contB[0]]=lar; curidmapB[extB[0]]=lar-1;
					curidmapC[0]=lar; curidmapC[2]=4-lar; curidmapC[1]=4-lar-1; curidmapC[3]=lar-1;	 
				}
				else
				{
					curidmapA[extA[2]]=4-lar; curidmapA[extA[1]]=lar-1; curidmapA[contA[0]]=4-lar-1;
					curidmapB[contB[0]]=4-lar; curidmapB[extB[0]]=4-lar-1;
					curidmapC[0]=lar; curidmapC[2]=4-lar; curidmapC[1]=lar-1; curidmapC[3]=4-lar-1;
				}
			}
		}
		else if(A->tensor_str.at(extA[0])  ==  A->tensor_str.at(extA[1]))
		{
			curidmapA[extA[1]]=lar-1; curidmapA[extA[2]]=4-lar; curidmapA[contA[0]]=4-lar-1;
			curidmapB[contB[0]]=4-lar; curidmapB[extB[0]]=4-lar-1;
			curidmapC[0]=lar; curidmapC[1]=lar-1; curidmapC[2]=4-lar; curidmapC[3]=4-lar-1;
		}
		else if(A->tensor_str.at(extA[0])  ==  A->tensor_str.at(extA[2]))
		{
			curidmapA[extA[2]]=lar-1; curidmapA[extA[1]]=4-lar; curidmapA[contA[0]]=4-lar-1;
			curidmapB[contB[0]]=4-lar; curidmapB[extB[0]]=4-lar-1;
			curidmapC[0]=lar; curidmapC[2]=lar-1; curidmapC[1]=4-lar; curidmapC[3]=4-lar-1;
		}
		else
		{
			curidmapA[contA[0]]=lar-1; curidmapA[extA[1]]=4-lar; curidmapA[extA[2]]=4-lar-1;
			curidmapB[contB[0]]=lar; curidmapB[extB[0]]=lar-1;
			curidmapC[0]=lar; curidmapC[2]=4-lar; curidmapC[1]=4-lar-1; curidmapC[3]=lar-1;
		}	
	}

}	





void cost::contract_4D_4D_2D(int* &idmapA, int* &idmapB, int* &idmapC, int* &ndimG, int* &npgrid)
{	
	
	ndimG[0]=4; //Always 4D grid for these contractions
	npgrid=new int[ndimG[0]];
	
	int fsize=2*ceil(sqrt(nprocs));
	int* temp = new int[fsize];
	int nfacts=0;

	for(int i=1; i<ceil(sqrt(nprocs));i++) if(nprocs%i==0) {temp[nfacts]=i; nfacts++; temp[nfacts]=nprocs/i; nfacts++;}
	int *facts = new int[4*nfacts];
	for(int i=0;i<nfacts;i++){for(int j=0; j<4; j++)facts[i*4+j]=temp[i];}
	
	int* curgrid = new int[4];
	found=false;
	Grid_Comb_4Dc(facts,curgrid,0,4*nfacts-1,0,4, idmapA, idmapB, idmapC, npgrid);
}


void cost::Grid_Comb_4Dc(int arr[], int data[], int start, int end, int index, int r,int* &idmapA, int* &idmapB, int* &idmapC, int* &npgrid)
{	
	if(found)return;
	if(index==r)
	{
		int temp=1;
		for (int i=0;i<4;i++) temp=temp*data[i];
		if(temp==nprocs && data[0]==data[1] && data[2]==data[3]) //Additional conditions to ensure symmetric dimensions are in same sized dimensions
		{	
			best_Cost_GGrid_4Dc(data, idmapA, idmapB, idmapC);	     	
			for(int i=0;i<4;i++) npgrid[i]=data[i];
			found=true;
		}
	return;
	}

	for (int i=start; i<=end && end-i+1 >= r-index; i++)
    	{
        	data[index] = arr[i];
        	Grid_Comb_4Dc(arr, data, i+1, end, index+1, r, idmapA, idmapB, idmapC, npgrid);
  
        	// Remove duplicates
        	while (arr[i] == arr[i+1])
             	i++;
    	}

}


void cost::best_Cost_GGrid_4Dc(int* &curGrid, int* &curidmapA, int* &curidmapB, int* &curidmapC)
{	
	
	curidmapA = new int[4];
	curidmapB=new int[4];
	curidmapC=new int[2];
	int contA[3]={};
	int extA[1]={};
	int contB[3]={};
	int extB[1]={};

	int ccnt=0, ecnt=0;
	for(int i=0; i<dimA; i++)
	    	{
        		for(int j=0; j<dimB; j++)
       		 	{
        	    		if(A->contr_dim_str[i].compare(B->contr_dim_str[j]) == 0)
          			{
               		 		// i in A and j in B are contracting indices
                			contA[ccnt]=i; contB[ccnt]=j;
					ccnt++;
						
            			}
        		}
    		}

	for(int i=0;i<4;i++) {if(i!=contA[0]) extA[ecnt]=i; ecnt++;}
	ecnt=0;
	for(int i=0;i<4;i++) {if(i!=contB[0]) extB[ecnt]=i; ecnt++;}

	int lar=1;
	if(curGrid[3] > curGrid[1]) lar=3;

	int size_A=1, size_B=1;
	for(int i=0;i<4;i++){ size_A*=A->tensor_size[i];	size_B*=B->tensor_size[i];}
	
	bool symA1f=false,symA2f=false,symB1f=false,symB2f=false;
	
	
	if(A->tensor_str.at(extA[0]) == 'a' || A->tensor_str.at(contA[0]) == 'A' || A->tensor_str.at(contA[1]) == 'A' || A->tensor_str.at(contA[2]) == 'A') symA1f=true;
	if(A->tensor_str.at(extA[0]) == 'b' || A->tensor_str.at(contA[0]) == 'B' || A->tensor_str.at(contA[1]) == 'B' || A->tensor_str.at(contA[2]) == 'B') symA2f=true;
	if(B->tensor_str.at(extA[0]) == 'a' || B->tensor_str.at(contA[0]) == 'A' || B->tensor_str.at(contA[1]) == 'A' || B->tensor_str.at(contA[2]) == 'A') symB1f=true;
	if(B->tensor_str.at(extA[0]) == 'b' || B->tensor_str.at(contA[0]) == 'B' || B->tensor_str.at(contA[1]) == 'B' || B->tensor_str.at(contA[2]) == 'B') symB2f=true;
	
	
	if(size_A >= size_B)
	{	
		curidmapA[extA[0]]=lar;
		if(!symA1f && !symA2f)
		{
			curidmapA[contA[0]]=lar-1; curidmapA[contA[1]]=4-lar; curidmapA[contA[2]]=4-lar-1;
		}
		else if(symA1f && A->tensor_str.at(extA[0]) == 'a')
		{
			if(A->tensor_str.at(contA[0])=='A'){ curidmapA[contA[0]]=lar-1; curidmapA[contA[1]]=4-lar; curidmapA[contA[2]]=4-lar-1; }
			else if(A->tensor_str.at(contA[1])=='A'){ curidmapA[contA[1]]=lar-1; curidmapA[contA[0]]=4-lar; curidmapA[contA[2]]=4-lar-1; }
			else if(A->tensor_str.at(contA[2])=='A'){ curidmapA[contA[2]]=lar-1; curidmapA[contA[0]]=4-lar; curidmapA[contA[0]]=4-lar-1; }
		}
		else if(symA2f && A->tensor_str.at(extA[0]) == 'b')
		{
			if(A->tensor_str.at(contA[0])=='B'){ curidmapA[contA[0]]=lar-1; curidmapA[contA[1]]=4-lar; curidmapA[contA[2]]=4-lar-1; }
			else if(A->tensor_str.at(contA[1])=='B'){ curidmapA[contA[1]]=lar-1; curidmapA[contA[0]]=4-lar; curidmapA[contA[2]]=4-lar-1; }
			else if(A->tensor_str.at(contA[2])=='B'){ curidmapA[contA[2]]=lar-1; curidmapA[contA[0]]=4-lar; curidmapA[contA[0]]=4-lar-1; }
		}
		else if(symA1f && A->tensor_str.at(contA[0]) == 'A')
		{
			if(A->tensor_str.at(contA[1]) == 'A'){curidmapA[contA[0]]=4-lar; curidmapA[contA[1]]=4-lar-1; curidmapA[contA[2]]=lar-1;}
			else if(A->tensor_str.at(contA[2]) == 'A'){curidmapA[contA[0]]=4-lar; curidmapA[contA[2]]=4-lar-1; curidmapA[contA[1]]=lar-1;}
		}
		else if(symA2f && A->tensor_str.at(contA[0]) == 'B')
		{
			if(A->tensor_str.at(contA[1]) == 'B'){curidmapA[contA[0]]=4-lar; curidmapA[contA[1]]=4-lar-1; curidmapA[contA[2]]=lar-1;}
			else if(A->tensor_str.at(contA[2]) == 'B'){curidmapA[contA[0]]=4-lar; curidmapA[contA[2]]=4-lar-1; curidmapA[contA[1]]=lar-1;}
		}
		else if(symA1f && A->tensor_str.at(contA[1]) == 'A')
		{
			curidmapA[contA[1]]=4-lar; curidmapA[contA[2]]=4-lar-1; curidmapA[contA[0]]=lar-1;
		}
		else if(symA2f && 	A->tensor_str.at(contA[1]) == 'B')
		{
			curidmapA[contA[1]]=4-lar; curidmapA[contA[2]]=4-lar-1; curidmapA[contA[0]]=lar-1;
		}
		
		if(!symB1f && !symB2f)
		{
			curidmapB[extB[0]]=lar-1; curidmapB[contB[0]]=lar; curidmapB[contB[1]]=4-lar; curidmapB[contB[2]]=4-lar-1;
			curidmapC[0]=lar; curidmapC[1]=lar-1;
		}
		else if(symB1f && B->tensor_str.at(extB[0]) == 'a')
		{
			curidmapB[extB[0]]=4-lar; 
			if(B->tensor_str.at(contB[0]) == 'A'){ curidmapB[contB[0]]=4-lar-1; curidmapB[contB[1]]=lar; curidmapB[contB[2]]=lar-1; curidmapC[0]=lar; curidmapC[1]=4-lar;}
			else if(B->tensor_str.at(contB[1]) == 'A'){ curidmapB[contB[1]]=4-lar-1; curidmapB[contB[0]]=lar; curidmapB[contB[2]]=lar-1; curidmapC[0]=lar; curidmapC[1]=4-lar;}
			else if(B->tensor_str.at(contB[2]) == 'A'){ curidmapB[contB[2]]=4-lar-1; curidmapB[contB[1]]=lar; curidmapB[contB[0]]=lar-1; curidmapC[0]=lar; curidmapC[1]=4-lar;}
		}	
		else if	(symB2f && B->tensor_str.at(extB[0]) == 'b')
		{
			curidmapB[extB[0]]=4-lar; 
			if(B->tensor_str.at(contB[0]) == 'B'){ curidmapB[contB[0]]=4-lar-1; curidmapB[contB[1]]=lar; curidmapB[contB[2]]=lar-1; curidmapC[0]=lar; curidmapC[1]=4-lar;}
			else if(B->tensor_str.at(contB[1]) == 'B'){ curidmapB[contB[1]]=4-lar-1; curidmapB[contB[0]]=lar; curidmapB[contB[2]]=lar-1; curidmapC[0]=lar; curidmapC[1]=4-lar;}
			else if(B->tensor_str.at(contB[2]) == 'B'){ curidmapB[contB[2]]=4-lar-1; curidmapB[contB[1]]=lar; curidmapB[contB[0]]=lar-1; curidmapC[0]=lar; curidmapC[1]=4-lar;}
		}
		else if	(symB1f && B->tensor_str.at(contB[0]) == 'A')
		{
			curidmapB[extB[0]]=lar-1; 
			if(B->tensor_str.at(contB[1]) == 'A'){ curidmapB[contB[0]]=4-lar; curidmapB[contB[1]]=4-lar-1; curidmapB[contB[2]]=lar-1; curidmapC[0]=lar; curidmapC[1]=lar-1;}
			else if(B->tensor_str.at(contB[2]) == 'A'){ curidmapB[contB[1]]=lar-1; curidmapB[contB[0]]=4-lar; curidmapB[contB[2]]=4-lar-1; curidmapC[0]=lar; curidmapC[1]=lar-1;}
		}	
		else if	(symB1f && B->tensor_str.at(contB[1]) == 'A')
		{
			curidmapB[extB[0]]=lar-1; 
			 curidmapB[contB[0]]=lar; curidmapB[contB[1]]=4-lar; curidmapB[contB[2]]=4-lar-1; curidmapC[0]=lar; curidmapC[1]=lar-1;
		}
		else if	(symB2f && B->tensor_str.at(contB[0]) == 'B')
		{
			curidmapB[extB[0]]=lar-1; 
			if(B->tensor_str.at(contB[1]) == 'B'){ curidmapB[contB[0]]=4-lar; curidmapB[contB[1]]=4-lar-1; curidmapB[contB[2]]=lar-1; curidmapC[0]=lar; curidmapC[1]=lar-1;}
			else if(B->tensor_str.at(contB[2]) == 'B'){ curidmapB[contB[1]]=lar-1; curidmapB[contB[0]]=4-lar; curidmapB[contB[2]]=4-lar-1; curidmapC[0]=lar; curidmapC[1]=lar-1;}
		}	
		else if	(symB2f && B->tensor_str.at(contB[1]) == 'B')
		{
			curidmapB[extB[0]]=lar-1; 
			 curidmapB[contB[0]]=lar; curidmapB[contB[1]]=4-lar; curidmapB[contB[2]]=4-lar-1; curidmapC[0]=lar; curidmapC[1]=lar-1;
		}
	}
	else
	{
		curidmapB[extB[0]]=lar;
		if(!symB1f && !symB2f)
		{
			curidmapB[contB[0]]=lar-1; curidmapB[contB[1]]=4-lar; curidmapB[contB[2]]=4-lar-1;
		}
		else if(symB1f && B->tensor_str.at(extB[0]) == 'a')
		{
			if(B->tensor_str.at(contB[0])=='A'){ curidmapB[contB[0]]=lar-1; curidmapB[contB[1]]=4-lar; curidmapB[contB[2]]=4-lar-1; }
			else if(B->tensor_str.at(contB[1])=='A'){ curidmapB[contB[1]]=lar-1; curidmapB[contB[0]]=4-lar; curidmapB[contB[2]]=4-lar-1; }
			else if(B->tensor_str.at(contB[2])=='A'){ curidmapB[contB[2]]=lar-1; curidmapB[contB[0]]=4-lar; curidmapB[contB[0]]=4-lar-1; }
		}
		else if(symB2f && B->tensor_str.at(extB[0]) == 'b')
		{
			if(B->tensor_str.at(contB[0])=='B'){ curidmapB[contB[0]]=lar-1; curidmapB[contB[1]]=4-lar; curidmapB[contB[2]]=4-lar-1; }
			else if(B->tensor_str.at(contB[1])=='B'){ curidmapB[contB[1]]=lar-1; curidmapB[contB[0]]=4-lar; curidmapB[contB[2]]=4-lar-1; }
			else if(B->tensor_str.at(contB[2])=='B'){ curidmapB[contB[2]]=lar-1; curidmapB[contB[0]]=4-lar; curidmapB[contB[0]]=4-lar-1; }
		}
		else if(symB1f && B->tensor_str.at(contB[0]) == 'A')
		{
			if(B->tensor_str.at(contB[1]) == 'A'){curidmapB[contB[0]]=4-lar; curidmapB[contB[1]]=4-lar-1; curidmapB[contB[2]]=lar-1;}
			else if(B->tensor_str.at(contB[2]) == 'A'){curidmapB[contB[0]]=4-lar; curidmapB[contB[2]]=4-lar-1; curidmapB[contB[1]]=lar-1;}
		}
		else if(symB2f && B->tensor_str.at(contB[0]) == 'B')
		{
			if(B->tensor_str.at(contB[1]) == 'B'){curidmapB[contB[0]]=4-lar; curidmapB[contB[1]]=4-lar-1; curidmapB[contB[2]]=lar-1;}
			else if(B->tensor_str.at(contB[2]) == 'B'){curidmapB[contB[0]]=4-lar; curidmapB[contB[2]]=4-lar-1; curidmapB[contB[1]]=lar-1;}
		}
		else if(symB1f && B->tensor_str.at(contB[1]) == 'A')
		{
			curidmapB[contB[1]]=4-lar; curidmapB[contB[2]]=4-lar-1; curidmapB[contB[0]]=lar-1;
		}
		else if(symB2f && 	B->tensor_str.at(contB[1]) == 'B')
		{
			curidmapB[contB[1]]=4-lar; curidmapB[contB[2]]=4-lar-1; curidmapB[contB[0]]=lar-1;
		}
		
		if(!symA1f && !symA2f)
		{
			curidmapA[extA[0]]=lar-1; curidmapA[contA[0]]=lar; curidmapA[contA[1]]=4-lar; curidmapA[contA[2]]=4-lar-1;
			curidmapC[0]=lar; curidmapC[1]=lar-1;
		}
		else if(symA1f && A->tensor_str.at(extA[0]) == 'a')
		{
			curidmapA[extA[0]]=4-lar; 
			if(A->tensor_str.at(contA[0]) == 'A'){ curidmapA[contA[0]]=4-lar-1; curidmapA[contA[1]]=lar; curidmapA[contA[2]]=lar-1; curidmapC[0]=lar; curidmapC[1]=4-lar;}
			else if(A->tensor_str.at(contA[1]) == 'A'){ curidmapA[contA[1]]=4-lar-1; curidmapA[contA[0]]=lar; curidmapA[contA[2]]=lar-1; curidmapC[0]=lar; curidmapC[1]=4-lar;}
			else if(A->tensor_str.at(contA[2]) == 'A'){ curidmapA[contA[2]]=4-lar-1; curidmapA[contA[1]]=lar; curidmapA[contA[0]]=lar-1; curidmapC[0]=lar; curidmapC[1]=4-lar;}
		}	
		else if	(symA2f && A->tensor_str.at(extA[0]) == 'b')
		{
			curidmapA[extA[0]]=4-lar; 
			if(A->tensor_str.at(contA[0]) == 'B'){ curidmapA[contA[0]]=4-lar-1; curidmapA[contA[1]]=lar; curidmapA[contA[2]]=lar-1; curidmapC[0]=lar; curidmapC[1]=4-lar;}
			else if(A->tensor_str.at(contA[1]) == 'B'){ curidmapA[contA[1]]=4-lar-1; curidmapA[contA[0]]=lar; curidmapA[contA[2]]=lar-1; curidmapC[0]=lar; curidmapC[1]=4-lar;}
			else if(A->tensor_str.at(contA[2]) == 'B'){ curidmapA[contA[2]]=4-lar-1; curidmapA[contA[1]]=lar; curidmapA[contA[0]]=lar-1; curidmapC[0]=lar; curidmapC[1]=4-lar;}
		}
		else if	(symA1f && A->tensor_str.at(contA[0]) == 'A')
		{
			curidmapA[extA[0]]=lar-1; 
			if(A->tensor_str.at(contA[1]) == 'A'){ curidmapA[contA[0]]=4-lar; curidmapA[contA[1]]=4-lar-1; curidmapA[contA[2]]=lar-1; curidmapC[0]=lar; curidmapC[1]=lar-1;}
			else if(A->tensor_str.at(contA[2]) == 'A'){ curidmapA[contA[1]]=lar-1; curidmapA[contA[0]]=4-lar; curidmapA[contA[2]]=4-lar-1; curidmapC[0]=lar; curidmapC[1]=lar-1;}
		}	
		else if	(symA1f && A->tensor_str.at(contA[1]) == 'A')
		{
			curidmapA[extA[0]]=lar-1; 
			 curidmapA[contA[0]]=lar; curidmapA[contA[1]]=4-lar; curidmapA[contA[2]]=4-lar-1; curidmapC[0]=lar; curidmapC[1]=lar-1;
		}
		else if	(symA2f && A->tensor_str.at(contA[0]) == 'B')
		{
			curidmapA[extA[0]]=lar-1; 
			if(A->tensor_str.at(contA[1]) == 'B'){ curidmapA[contA[0]]=4-lar; curidmapA[contA[1]]=4-lar-1; curidmapA[contA[2]]=lar-1; curidmapC[0]=lar; curidmapC[1]=lar-1;}
			else if(A->tensor_str.at(contA[2]) == 'B'){ curidmapA[contA[1]]=lar-1; curidmapA[contA[0]]=4-lar; curidmapA[contA[2]]=4-lar-1; curidmapC[0]=lar; curidmapC[1]=lar-1;}
		}	
		else if	(symA2f && A->tensor_str.at(contA[1]) == 'B')
		{
			curidmapA[extA[0]]=lar-1; 
			 curidmapA[contA[0]]=lar; curidmapA[contA[1]]=4-lar; curidmapA[contA[2]]=4-lar-1; curidmapC[0]=lar; curidmapC[1]=lar-1;
		}
	}	
}	



void cost::contract_2D_4D_4D(int* &idmapA, int* &idmapB, int* &idmapC, int* &ndimG, int* &npgrid)
{	
	
	ndimG[0]=4; //Always 4D grid for these contractions
	npgrid=new int[ndimG[0]];
	
	int fsize=2*ceil(sqrt(nprocs));
	int* temp = new int[fsize];
	int nfacts=0;

	for(int i=1; i<ceil(sqrt(nprocs));i++) if(nprocs%i==0) {temp[nfacts]=i; nfacts++; temp[nfacts]=nprocs/i; nfacts++;}
	int *facts = new int[4*nfacts];
	for(int i=0;i<nfacts;i++){for(int j=0; j<4; j++)facts[i*4+j]=temp[i];}
	
	int* curgrid = new int[4];
	found=false;
	Grid_Comb_4Dd(facts,curgrid,0,4*nfacts-1,0,4, idmapA, idmapB, idmapC, npgrid);
}



void cost::Grid_Comb_4Dd(int arr[], int data[], int start, int end, int index, int r,int* &idmapA, int* &idmapB, int* &idmapC, int* &npgrid)
{	
	if(found)return;
	if(index==r)
	{
		int temp=1;
		for (int i=0;i<4;i++) temp=temp*data[i];
		if(temp==nprocs && data[0]==data[1] && data[2]==data[3]) //Additional conditions to ensure symmetric dimensions are in same sized dimensions
		{	
			best_Cost_GGrid_4Dd(data, idmapA, idmapB, idmapC);	     	
			for(int i=0;i<4;i++) npgrid[i]=data[i];
			found=true;
		}
	return;
	}

	for (int i=start; i<=end && end-i+1 >= r-index; i++)
    	{
        	data[index] = arr[i];
        	Grid_Comb_4Dd(arr, data, i+1, end, index+1, r, idmapA, idmapB, idmapC, npgrid);
  
        	// Remove duplicates
        	while (arr[i] == arr[i+1])
             	i++;
    	}

}


void cost::best_Cost_GGrid_4Dd(int* &curGrid, int* &curidmapA, int* &curidmapB, int* &curidmapC)
{	
	
	curidmapA = new int[2];
	curidmapB=new int[4];
	curidmapC=new int[4];
	int contB[1]={};
	int extB[3]={};
	int contA[1]={};
	int extA[1]={};

	int ccnt=0, ecnt=0;
	for(int i=0; i<dimA; i++)
	    	{
        		for(int j=0; j<dimB; j++)
       		 	{
        	    		if(A->contr_dim_str[i].compare(B->contr_dim_str[j]) == 0)
          			{
               		 		// i in A and j in B are contracting indices
                			contA[ccnt]=i; contB[ccnt]=j;
					ccnt++;
						
            			}
        		}
    		}

	for(int i=0;i<4;i++) {if(i!=contA[0]) extA[ecnt]=i; ecnt++;}
	ecnt=0;
	for(int i=0;i<4;i++) {if(i!=contB[0]) extB[ecnt]=i; ecnt++;}

	int lar=1;
	if(curGrid[3] > curGrid[1]) lar=3;

	int size_A=1, size_B=1;
	for(int i=0;i<4;i++){ size_A*=A->tensor_size[i];	size_B*=B->tensor_size[i];}
	
	int temp;
	for(int i=1; i<3; i++){if(B->tensor_size[extB[i]] > B->tensor_size[extB[0]] ) {temp = extB[0]; extA[0]=extB[i]; extB[i]=temp;}}
	
	
	if(B->tensor_size[extB[0]] >= B->tensor_size[extB[1]] && B->tensor_size[extB[0]] >= B->tensor_size[extB[2]])
	{
		curidmapB[extB[0]]=lar;
		if(B->tensor_str.at(extB[0]) == 'c' || B->tensor_str.at(extB[0]) == 'd')
		{ 
			if(B->tensor_size[extB[1]] >= B->tensor_size[extB[2]])
			{
				if(B->tensor_str.at(extB[1]) == 'c' || B->tensor_str.at(extB[1]) == 'd')
				{
					curidmapB[extB[1]]=lar-1; curidmapB[extB[2]]=4-lar; curidmapB[contB[0]]=4-lar-1;
					curidmapA[contA[0]]=4-lar; curidmapA[extA[0]]=4-lar-1;
					curidmapC[0]=lar; curidmapC[1]=lar-1; curidmapC[2]=4-lar; curidmapC[3]=4-lar-1;
				}
				else if(B->tensor_str.at(extB[1]) == B->tensor_str.at(extB[2]))
				{
					curidmapB[extB[1]]=4-lar; curidmapB[extB[2]]=4-lar-1; curidmapB[contB[0]]=lar-1;
					curidmapA[contA[0]]=lar; curidmapA[extA[0]]=lar-1;
					curidmapC[0]=lar; curidmapC[1]=4-lar; curidmapC[2]=4-lar-1; curidmapC[3]=lar-1;	 
				}
				else
				{
					curidmapB[extB[1]]=4-lar; curidmapB[extB[2]]=lar-1; curidmapB[contB[0]]=4-lar-1;
					curidmapA[contA[0]]=4-lar; curidmapA[extA[0]]=4-lar-1;
					curidmapC[0]=lar; curidmapC[1]=4-lar; curidmapC[2]=lar-1; curidmapC[3]=4-lar-1;
				}
			}
			else if(B->tensor_size[extB[2]] >= B->tensor_size[extB[1]])
			{
				if(B->tensor_str.at(extB[2]) == 'c' || B->tensor_str.at(extB[2]) == 'd')
				{
					curidmapB[extB[2]]=lar-1; curidmapB[extB[1]]=4-lar; curidmapB[contB[0]]=4-lar-1;
					curidmapA[contA[0]]=4-lar; curidmapA[extA[0]]=4-lar-1;
					curidmapC[0]=lar; curidmapC[2]=lar-1; curidmapC[1]=4-lar; curidmapC[3]=4-lar-1;
				}
				else if(A->tensor_str.at(extA[2]) == A->tensor_str.at(extA[1]))
				{
					curidmapB[extB[2]]=4-lar; curidmapB[extB[1]]=4-lar-1; curidmapB[contB[0]]=lar-1;
					curidmapA[contA[0]]=lar; curidmapA[extA[0]]=lar-1;
					curidmapC[0]=lar; curidmapC[2]=4-lar; curidmapC[1]=4-lar-1; curidmapC[3]=lar-1;	 
				}
				else
				{
					curidmapB[extB[2]]=4-lar; curidmapB[extB[1]]=lar-1; curidmapB[contB[0]]=4-lar-1;
					curidmapA[contA[0]]=4-lar; curidmapA[extA[0]]=4-lar-1;
					curidmapC[0]=lar; curidmapC[2]=4-lar; curidmapC[1]=lar-1; curidmapC[3]=4-lar-1;
				}
			}
		}
		else if(B->tensor_str.at(extB[0])  ==  B->tensor_str.at(extB[1]))
		{
			curidmapB[extB[1]]=lar-1; curidmapB[extB[2]]=4-lar; curidmapB[contB[0]]=4-lar-1;
			curidmapA[contA[0]]=4-lar; curidmapA[extA[0]]=4-lar-1;
			curidmapC[0]=lar; curidmapC[1]=lar-1; curidmapC[2]=4-lar; curidmapC[3]=4-lar-1;
		}
		else if(B->tensor_str.at(extB[0])  ==  B->tensor_str.at(extB[2]))
		{
			curidmapB[extB[2]]=lar-1; curidmapB[extB[1]]=4-lar; curidmapB[contB[0]]=4-lar-1;
			curidmapA[contA[0]]=4-lar; curidmapA[extA[0]]=4-lar-1;
			curidmapC[0]=lar; curidmapC[2]=lar-1; curidmapC[1]=4-lar; curidmapC[3]=4-lar-1;
		}
		else
		{
			curidmapB[contB[0]]=lar-1; curidmapB[extB[1]]=4-lar; curidmapB[extB[2]]=4-lar-1;
			curidmapA[contA[0]]=lar; curidmapA[extA[0]]=lar-1;
			curidmapC[0]=lar; curidmapC[2]=4-lar; curidmapC[1]=4-lar-1; curidmapC[3]=lar-1;
		}	
	}

}	

/*
double cost::calcCost(int dimsM, int dimsB, int dimsC, Tensor* &A, Tensor* &B, Tensor* &C, int griddims, Grid* &G)
{
    
    list<pair<int,int>> contr_list = list<pair<int,int>>();
    list<pair<int,int>> DDO_list = list<pair<int,int>>();
    list<pair<int,int>> DDA_list = list<pair<int,int>>();
    int* ddopairs = new int[griddims];


//Find the list of contraction indices pairs
    for(int i=0; i<dimsA; i++)
    {
        for(int j=0; j<dimsB; j++)
        {
            if(A->contr_dim_str[i].compare(B->contr_dim_str[j]) == 0)
            {
                // i in A and j in B are contracting indices
                pair<int,int> p(i, j);
                contr_list.push_back(p);
                A->cntr_map[i] = 1;
                B->cntr_map[j] = 1;
            }
        }
    }


//initialize array that will hold the dims that will be reduced, broadcast and rotated.
    int* reduction_dims = new int[griddims];
    int* broadcast_dims = new int[griddims];
    int* rotate_dims = new int[griddims];


    memset(reduction_dims, 0, griddims*sizeof(int));
    memset(broadcast_dims, 0, griddims*sizeof(int));
    memset(rotate_dims, 0, griddims*sizeof(int));
   // memset(ddopairs, 99, griddims*sizeof(int));
    int* bindices= new int[dimsB];
    memset(bindices, 0, dimsB*sizeof(int));

for(int i=0; i<griddims;i++) ddopairs[i]=9999;

    for (std::list<pair<int,int>>::iterator cntr_it=contr_list.begin(); cntr_it != contr_list.end(); ++cntr_it)
    {
        int da = A->index_dimension_map[(*cntr_it).first];
        int db = B->index_dimension_map[(*cntr_it).second] ;
        bindices[(*cntr_it).second]=1;
        if(da == db)
        {
            DDA_list.push_back(*cntr_it);
            reduction_dims[da] = 1;
        }
        else
        {
            DDO_list.push_back(*cntr_it);
            broadcast_dims[da] = 1;
            broadcast_dims[db] = 1;
            ddopairs[da]=db;
        }
    }

    


    
     for(int i=0; i<dimsA; i++)
    {
        for(int j=0; j<dimsB; j++)
        {
            if(A->contr_dim_str[i].compare(B->contr_dim_str[j]) != 0 && 
           !A->cntr_map[i] && !B->cntr_map[j])
            {

        if(A->index_dimension_map[i] == B->index_dimension_map[j] && A->index_dimension_map[i] < griddims)
        {
            
            
            rotate_dims[A->index_dimension_map[i]]=1;
           
        }
            }
        }
    }

    
    int ts=0,tw=1;        
    int num_rotation=1;
    for(int i=0; i<griddims; i++)
    {
        if(rotate_dims[i]) num_rotation= num_rotation*G->pgrid[i];
    }


    
    double pNA=1,pNB=1,pNC=1;
    double SA=1,SB=1,SC=1;
    double rec_cost=0;

    for(int i=0; i<dimsA; i++)
    {
        if(A->index_dimension_map[i]<griddims) pNA = pNA * G->pgrid[A->index_dimension_map[i]];
        SA = SA * A->tensor_size[i];
    }
           


    for(int i=0; i<dimsB; i++)
    {
        if(B->index_dimension_map[i]<griddims) pNB = pNB * G->pgrid[B->index_dimension_map[i]];
        SB = SB * B->tensor_size[i];
    }
           

    for(int i=0; i<dimsC; i++)
    {
        if(C->index_dimension_map[i]<griddims) pNC = pNC * G->pgrid[C->index_dimension_map[i]];
        SC = SC * C->tensor_size[i];
    }
           
    double mA=SA/pNA;
    double mB=SB/pNB;
    double mC=SC/pNC;
    double gammam=mC;


    //computation product of all external indicies and contraction divided by the no. of procs across which they are distributed times 2

    mA = mA*sizeof(double);
    mB = mB*sizeof(double);
    mC = mC*sizeof(double);

    int red_procs=1;
    for(int i=0; i<griddims; i++)
    {
        if(reduction_dims[i]) red_procs = red_procs*G->pgrid[i];
    }

    
double br_cost=0;

    br_cost = bcost(0,mA,mB,tw,G,griddims,ddopairs,'A') + bcost(0,mA,mB,tw,G,griddims,ddopairs,'B');
   

     double red_cost = log2(red_procs)*(ts+tw*mC+gammam);

    double comp_cost=2;

    for(int i=0; i<dimsA; i++)
        {
            comp_cost=comp_cost*A->tensor_size[i];
            if(A->index_dimension_map[i]<griddims) comp_cost=comp_cost/G->pgrid[A->index_dimension_map[i]];
        }

    for(int i=0; i<dimsB; i++)
        {
            if(!bindices[i]) comp_cost=comp_cost*B->tensor_size[i];
            if(B->index_dimension_map[i]<grid_dims) comp_cost=comp_cost/G->pgrid[B->index_dimension_map[i]];
        }    

     double tot_cost = num_rotation * (br_cost + comp_cost + red_cost + ts + (tw*mA));

     return tot_cost;

}


bool cost::check_redistr(Tensor* &T, Tensor* &C, string* &t, string* &c)
{

    // Initialize new idmap to the default value (old idmap)
    int* new_idmap = new int[T->dims];
    memcpy(new_idmap, T->index_dimension_map, T->dims * sizeof(int)); 

    bool redistr_flag = false;

    for(int i=0; i < C->dims; i++)
    {
        for(int k =0; k < T->dims; k++)
        {
            // Find the external dimension k in T with the analogous dimension i in C
            if(C->contr_dim_str[i].compare(T->contr_dim_str[k]) == 0)
            {
                // Their index maps should be same... 
                // If not, redistribute T so that both of their physical dimension is same
                if(C->index_dimension_map[i] != new_idmap[k]  &&  new_idmap[k] < grid_dims)
                {
                    redistr_flag = true;

                    //if(rank==0) cout<< C->contr_dim_str[i] << "  " << T->contr_dim_str[k] << " i= " << i << " k= "<<k << " c_idmap= " <<C->index_dimension_map[i] << " t_idmap= " <<T->index_dimension_map[k] << "  new_idmap[i] = " << new_idmap[i] << "  new_idmap[k]= " << new_idmap[k] <<endl;
                    break;
                }
            }
        }
    }
    //if(rank==0){ cout << endl << " new_idmap = " ; print_tile_addr(A->dims, new_idmap); cout << endl;}
    return redistr_flag; 
}


double cost::bcost(int iter,double MA,double MB,double tw, Grid* &G, int griddims, int* &ddopairs, char tnsr)
{   
    if(tnsr == 'A')
    for(int i=iter; i<griddims; i++)
        if(ddopairs[i]<999)
        {   
            return ((G->pgrid[i]*MA*2) + (G->pgrid[i]*bcost(i+1, MA, MB, tw, G, grid_dims, ddopairs, 'A')));
        }

        if(tnsr == 'B') 
    for(int i=iter; i<griddims; i++)
        if(ddopairs[i]<999)
        {
             return ((G->pgrid[ddopairs[i]]*MB*2) + (G1->pgrid[ddopairs[i]]*bcost(i+1, MA, MB, tw, G, grid_dims, ddopairs, 'B')));
        }

          
      return 0;    
}


double cost::getclock() {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return (tp.tv_sec + tp.tv_usec*1.0e-6);
}


void cost::best_Cost(int* &new_idmap_A,int* &new_idmap_B,int* &new_idmap_C,int* &new_griddims,int* &new_gridcount)
{


int* origg = new int[grid_dims];
for(int i=0; i<grid_dims; i++) origg[i]=G1->pgrid[i];

int nprocs=1;
	for(int i=0; i<	grid_dims; i++) nprocs = nprocs*G1->pgrid[i];

    int nfacts=0;
int fsize=2*ceil(sqrt(nprocs));
    int* temp = new int[fsize];
int maxpdim=max(max(dimA,dimB),dimC);
for(int i=1; i<ceil(sqrt(nprocs));i++) if(nprocs%i==0) {temp[nfacts]=i; nfacts++; temp[nfacts]=nprocs/i; nfacts++;}
int *fact = new int[grid_dims*nfacts];
for(int i=0;i<nfacts;i++){for(int j=0; j<maxpdim; j++)fact[i*maxpdim+j]=temp[i];}
std::sort(fact, fact + (grid_dims*nfacts));
int* curr_comb = new int[grid_dims];
changed=false;
for(int k=2;k<=maxpdim;k++)
{
Grid_Comb(fact, curr_comb, 0, (grid_dims*nfacts-1), 0, k,nprocs);
}
if(changed) for(int i=0; i<grid_dims; i++){  G1->pgrid[i]=G2->pgrid[i];}



}



void cost::Grid_Comb(int arr[], int data[], int start, int end, int index, int r, int nprocs)
{
	if(index==r)
	{
		int temp=1;
		for (int i=0;i<grid_dims;i++) temp=temp*data[i];
		if(temp==nprocs)
		{
			
				changed=false;
                best_Cost_GGrid();
				if(changed){for(int i=0; i<grid_dims; i++){  G2->pgrid[i]=G1->pgrid[i];}}
            
			
		}
		return;
	}

for (int i=start; i<=end && end-i+1 >= r-index; i++)
    {
        data[index] = arr[i];
        Grid_Comb(arr, data, i+1, end, index+1, r, nprocs);
 
 
        // Remove duplicates
        while (arr[i] == arr[i+1])
             i++;
    }

}



void cost::best_Cost_GGrid()
{

    int mgrid = grid_dims+max(dimA,max(dimB,grid_dims));
    int* gridA = new int[mgrid];
    int* gridB = new int[mgrid];
    int* gridC = new int[mgrid];

	int amin[dimA],bmin[dimB],cmin[dimC];
    bool flag=false;
	double minCost=999999;
	double cCost;




for(int i=0;i<dimA;i++) cout<<A1->contr_dim_str[i]<<' ';
    cout<<'\n';
for(int i=0;i<dimB;i++) cout<<B1->contr_dim_str[i]<<' ';
    cout<<'\n';
    for(int i=0;i<dimC;i++) cout<<C1->contr_dim_str[i]<<' ';
        cout<<'\n';

    for(int i=0;i<mgrid;i++)
    {
        if(i<grid_dims)
        {
            gridA[i]=i;
            gridB[i]=i;
            gridC[i]=i;  
        }  
        else 
        {
            gridA[i]=grid_dims;
            gridB[i]=grid_dims;
            gridC[i]=grid_dims;
        }
    }
int cnt=0;
    do
    {
        for(int Ai=0; Ai<dimA; Ai++)
        {
            A1->index_dimension_map[Ai]=gridA[Ai];
        }
           // if(check_space(A1,dimA,grid_dims,G1))
        if(true)
            {
                do
                {
                   for(int Bi=0; Bi<dimB; Bi++)
                   {
                    B1->index_dimension_map[Bi]=gridB[Bi];
                   }
                    //if(check_space(B1,dimB,grid_dims,G1))
                   if(true)
                    {
                        do
                        {

                            for(int Ci=0; Ci<dimC; Ci++)
                             {
                                C1->index_dimension_map[Ci]=gridC[Ci];
                              }
                        
                      //  if(check_space(C1,dimC,grid_dims,G1) &&  !check_redistr(A1, C1, A1->contr_dim_str, C1->contr_dim_str) 
                         //   && !check_redistr(B1, C1, B1->contr_dim_str, C1->contr_dim_str))
                              if(check_space(A1,dimC,grid_dims,G1))
                        {
                           if(!check_redistr(A1, C1, A1->contr_dim_str, C1->contr_dim_str)){
                            if(!check_redistr(B1, C1, B1->contr_dim_str, C1->contr_dim_str)){
                                    flag=true;
				cCost=calcCost(dimA, dimB, dimC, A1, B1, C1, grid_dims, G1);
				if(cCost<minCost){
									changed=true;
									for(int Ai=0; Ai<dimA; Ai++) amin[Ai] = A1->index_dimension_map[Ai];
									for(int Bi=0; Bi<dimB; Bi++) bmin[Bi] = B1->index_dimension_map[Bi];        
									for(int Ci=0; Ci<dimC; Ci++) cmin[Ci] = C1->index_dimension_map[Ci];}
                        }}}
                        }while(next_partial_permutation(gridC,(gridC+mgrid),dimC));
                   }
                    
                }while(next_partial_permutation(gridB,(gridB+mgrid),dimB));
            }
        
    }while( next_partial_permutation(gridA,(gridA+mgrid),dimA));


if(changed)
{
    memcpy(new_idmap_A,amin[Ai],dimsA*sizeof(int)) ;
    memcpy(new_idmap_B,bmin[Bi],dimsB*sizeof(int)) ;
    memcpy(new_idmap_C,cmin[Ci],dimsC*sizeof(int)) ;
}
f= flag;
}

bool cost::check_space(Tensor* &A, int dimsC, int griddims, Grid* &G)
{	
	double PROCMEM1 = 8*pow(2,30);
    double SA=1,pNA=1,SB=1,SC=1,pNB=1,pNC=1;
     for(int i=0; i<dimA; i++)
    {
        if(A1->index_dimension_map[i]<grid_dims) pNA = pNA * G->pgrid[A1->index_dimension_map[i]];
        SA = SA * A->tensor_size[i];
    }

    for(int i=0; i<dimB; i++)
    {
        if(B1->index_dimension_map[i]<grid_dims) pNB = pNB * G->pgrid[B1->index_dimension_map[i]];
        SB = SB * B1->tensor_size[i];
    }

     for(int i=0; i<dimC; i++)
    {
        if(C1->index_dimension_map[i]<grid_dims) pNC = pNC * G->pgrid[C1->index_dimension_map[i]];
        SC = SC * C1->tensor_size[i];
    }

    if(((SA/pNA + SB/pNB + SC/pNC)*sizeof(double)) < PROCMEM1) return true;
    return false;

}
*/


