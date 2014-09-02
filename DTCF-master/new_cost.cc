#include "new_cost.h"


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
	minCost=99999;
    }


void cost::bestCost(int* &idmapA, int* &idmapB, int* &idmapC, int* &ndimG, int* &npgrid)
{
	if(dimA==2 && dimB==2 && dimC==2)
		contract_2D_2D_2D(idmapA, idmapB, idmapC, ndimG, npgrid);

}

void cost::contract_2D_2D_2D(int* &idmapA, int* &idmapB, int* &idmapC, int* &ndimG, int* &npgrid)
{
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
	
	int curgrid = new int[4];
	
	Grid_Comb_4D(facts,curgrid,0,4*nfacts-1,0,4, idmapA, idmapB, idmapC, npgrid);
}


void cost::Grid_Comb_4D(int arr[], int data[], int start, int end, int index, int r,int* &idmapA, int* &idmapB, int* &idmapC, int* &npgrid)
{
	if(index==r)
	{
		int temp=1;
		for (int i=0;i<grid_dims;i++) temp=temp*data[i];
		if(temp==nprocs && data[0]==data[1] && data[2]==data[3]) //Additional conditions to ensure symmetric dimensions are in same sized dimensions
		{	
			int* curidmapA = new int[dimA]; int* curidmapB = new int[dimB]; int* curidmapC = new int[dimC];
			int curCost=best_Cost_GGrid_4D(data, curidmapA, curidmapB, curidmapC);	     	
			if(curCost<minCost)
			{ 
				minCost=curCost;	
				for(int i=0;i<4; i++) { npgrid[i]=data[i]; idmapA[i]=curidmapA[i]; idmapB[i]=curidmapB[i]; idmapC[i]=curidmapC[i]; }
			}
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


best_Cost_GGrid_4D(int* &curGrid, int* &curidmapA, int* &curidmapB, int* &curidmapC)
{	
	int contr_AB[4]={999,999,999,999};
	int contr_A[4]={};
	int contr_B[4]={};
	for(int i=0; i<dimA; i++)
	    	{
        		for(int j=0; j<dimB; j++)
       		 	{
        	    		if(A->contr_dim_str[i].compare(B->contr_dim_str[j]) == 0)
          			{
               		 		// i in A and j in B are contracting indices
                			contr_AB[i]=j; contr_A[i]=1; contr_B[j]=1;
						
            			}
        		}
    		}

int lar=0;
for(i=1;i<4;i++) {if (curGrid[i]>curGrid[i-1]) lar=i};

int size_A=1, size_B=1;
for(i=0;i<4;i++){ size_A*=A->size[i];	size_B*=B->size[i]}

if(size_A > size_B)
{	
	//First case
	if(contr_A[0]=0 && ((A->size[0] >= A->size[1] && contr_A[1]=0) || (A->size[0]>=A->size[2] && contr_A[2]=0)  || (A->size[0]>=A->size[3] && contr_A[3]=0))
	{
		curidmapA[0]=lar;
		if(A->SG_index_map[0]=A->SG_index_map[1]){ curidmapA[1]=lar-1 curidmapA[2]=4-lar; curidmapA[3]=4-lar-1;}
		else if(A->SG_index_map[0]=A->SG_index_map[2]){ curidmapA[2]=lar-1; curidmapA[1]=4-lar; curidmapA[3]=4-lar-1;}
		else if(A->SG_index_map[0]=A->SG_index_map[3]){ curidmapA[3]=lar-1; curidmapA[1]=4-lar; curidmapA[2]=4-lar-1;}	
		else if(contr_A[1]=0)
		{ 
			if(A->SG_index_map[1]!=SG_index_map[2] && A->SG_index_map[1]!=SG_index_map[3]) { curidmapA[1]=lar-1; curidmapA[2]=4-lar; curidmapA[3]=4-lar-1;}
			else if(A->SG_index_map[1]= SG_index_map[2]) {curidmapA[1]=4-lar; curidmapA[2]=4-lar-1; curidmapA[3]=lar-1;}
			else {curidmapA[1]=4-lar; curidmapA[3]=4-lar-1; curidmap[2]=lar-1;}
		}
		else if(contr_A[2]=0)
		{
		 	if(A->SG_index_map[2]!=SG_index_map[1] && A->SG_index_map[2]!=SG_index_map[3]) { curidmapA[2]=lar-1; curidmapA[3]=4-lar; curidmapA[4]=4-lar-1;}
			else if(A->SG_index_map[2]= SG_index_map[1]) {curidmapA[2]=4-lar; curidmapA[1]=4-lar-1; curidmapA[3]=lar-1;}
			else {curidmapA[2]=4-lar; curidmapA[3]=4-lar-1; curidmap[1]=lar-1;}
		}
		
		else if(contr_A[3]=0)
		{
		 	if(A->SG_index_map[3]!=SG_index_map[1] && A->SG_index_map[3]!=SG_index_map[2]) { curidmapA[3]=lar-1; curidmapA[1]=4-lar; curidmapA[2]=4-lar-1;}
			else if(A->SG_index_map[3]= SG_index_map[1]) {curidmapA[3]=4-lar; curidmapA[1]=4-lar-1; curidmapA[2]=lar-1;}
			else {curidmapA[3]=4-lar; curidmapA[2]=4-lar-1; curidmap[1]=lar-1;}
		}


	//Second case
	if(contr_A[0]=1 && ((A->size[1] >= A->size[0] && contr_A[0]=0) || (A->size[1]>=A->size[2] && contr_A[2]=0)  || (A->size[1]>=A->size[3] && contr_A[3]=0))
	{
		curidmapA[1]=lar;
		if(A->SG_index_map[0]=A->SG_index_map[1]){ curidmapA[1]=lar-1 curidmapA[2]=4-lar; curidmapA[3]=4-lar-1;}
		else if(A->SG_index_map[0]=A->SG_index_map[2]){ curidmapA[2]=lar-1; curidmapA[1]=4-lar; curidmapA[3]=4-lar-1;}
		else if(A->SG_index_map[0]=A->SG_index_map[3]){ curidmapA[3]=lar-1; curidmapA[1]=4-lar; curidmapA[2]=4-lar-1;}	
		else if(contr_A[1]=0)
		{ 
			if(A->SG_index_map[1]!=SG_index_map[2] && A->SG_index_map[1]!=SG_index_map[3]) { curidmapA[1]=lar-1; curidmapA[2]=4-lar; curidmapA[3]=4-lar-1;}
			else if(A->SG_index_map[1]= SG_index_map[2]) {curidmapA[1]=4-lar; curidmapA[2]=4-lar-1; curidmapA[3]=lar-1;}
			else {curidmapA[1]=4-lar; curidmapA[3]=4-lar-1; curidmap[2]=lar-1;}
		}
		else if(contr_A[2]=0)
		{
		 	if(A->SG_index_map[2]!=SG_index_map[1] && A->SG_index_map[2]!=SG_index_map[3]) { curidmapA[2]=lar-1; curidmapA[3]=4-lar; curidmapA[4]=4-lar-1;}
			else if(A->SG_index_map[2]= SG_index_map[1]) {curidmapA[2]=4-lar; curidmapA[1]=4-lar-1; curidmapA[3]=lar-1;}
			else {curidmapA[2]=4-lar; curidmapA[3]=4-lar-1; curidmap[1]=lar-1;}
		}
		
		else if(contr_A[3]=0)
		{
		 	if(A->SG_index_map[3]!=SG_index_map[1] && A->SG_index_map[3]!=SG_index_map[2]) { curidmapA[3]=lar-1; curidmapA[1]=4-lar; curidmapA[2]=4-lar-1;}
			else if(A->SG_index_map[3]= SG_index_map[1]) {curidmapA[3]=4-lar; curidmapA[1]=4-lar-1; curidmapA[2]=lar-1;}
			else {curidmapA[3]=4-lar; curidmapA[2]=4-lar-1; curidmap[1]=lar-1;}
		}
}	 



/*
double cost::calcCost(int dimsA, int dimsB, int dimsC, Tensor* &A, Tensor* &B, Tensor* &C, int griddims, Grid* &G)
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


