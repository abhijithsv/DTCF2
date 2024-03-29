#include "mpi.h"
#include "grid.h"
#include "tensor.h"
#include "contraction.h"
#include "redistribute.h"
#include "cost.h"
using namespace std;

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int rank,np;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    MPI_Comm_rank(MPI_COMM_WORLD,&np);

    int size = 128;
    int block_grid = 8;
    int grid_dim = 4;
    int* pgrid = new int[grid_dim];
    for(int i =0; i< grid_dim; i++)
    {
    pgrid[i] = 2;
    }
    if(rank==0) cout<<"Grid about to initialized"<<endl;

    Grid* grid = new Grid(grid_dim, pgrid);

    if(rank == 0) cout<<"Grid initialized"<<endl;

    int* size_A = new int[4];
    int* block_grid_A = new int[4];
    int* idmap_A = new int[4];

    size_A[0] = size;
    size_A[1] = size;
    size_A[2] = size;
    size_A[3] = size;

    block_grid_A[0] = block_grid;
    block_grid_A[1] = block_grid;
    block_grid_A[2] = block_grid;
    block_grid_A[3] = block_grid;

    idmap_A[0] = 0;
    idmap_A[1] = 1;
    idmap_A[2] = 2;
    idmap_A[3] = 3;

    Tensor* A = new Tensor("cccc", idmap_A, size_A, block_grid_A, grid);
    A->initialize();

    if(rank == 0) cout<<"Tensor A initialized"<<endl;

    int* size_B = new int[4];
    int* block_grid_B = new int[4];
    int* idmap_B = new int[4];

    size_B[0] = size;
    size_B[1] = size;
    size_B[2] = size;
    size_B[3] = size;

    block_grid_B[0] = block_grid;
    block_grid_B[1] = block_grid;
    block_grid_B[2] = block_grid;
    block_grid_B[3] = block_grid;

    idmap_B[0] = 2;
    idmap_B[1] = 3;
    idmap_B[2] = 4;
    idmap_B[3] = 4;

    Tensor* B = new Tensor("cccc", idmap_B, size_B, block_grid_B, grid);
    B->initialize();

    if(rank == 0) cout<<"Tensor B initialized"<<endl;

    int* size_C = new int[4];
    int* block_grid_C = new int[4];
    int* idmap_C = new int[4];

    size_C[0] = size;
    size_C[1] = size;
    size_C[2] = size;
    size_C[3] = size;

    block_grid_C[0] = block_grid;
    block_grid_C[1] = block_grid;
    block_grid_C[2] = block_grid;
    block_grid_C[3] = block_grid;

    idmap_C[0] = 0;
    idmap_C[1] = 1;
    idmap_C[2] = 2;
    idmap_C[3] = 3;

    Tensor* C = new Tensor("cccc", idmap_C, size_C, block_grid_C, grid);
    C->initialize();

    if(rank == 0) cout<<"Tensor C initialized"<<endl;

    double time =0.0;
   
    Contraction* C0 = new Contraction(A, B, C, grid);
    //time =-MPI_Wtime();
    C0->contract( "a,b,k,l", "m,n,k,l", "a,b,m,n");

    time +=MPI_Wtime();
    
    if(rank==0) cout << "Total Time ="<<time<<endl<<endl;
    MPI_Finalize();
}
