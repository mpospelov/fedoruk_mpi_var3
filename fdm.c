#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <assert.h>
#include <sys/time.h>
#include <mpi.h>

#define _REENTRANT

#define floatMPI MPI_DOUBLE
typedef double float_t;

void plotresult(float_t *A,float_t *Out, int N,int rank,int total, FILE* gnuplot);

int main(int argc, char **argv) {


  // typedef double float_t;
  int i,f,j,k;
  assert(argc==3);

  int N = atoi(argv[1]);

  int Tmax =  atoi(argv[2]);

  float_t *A;  //plate matrix
  float_t *Out; //just for root
  int n;    // Ширина ленты матрицы
    // Номер первой строки для вычисления

  int myrank, total;
  struct timeval startTime,endTime;

  gettimeofday(&startTime,NULL);

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &total);
  MPI_Comm_rank (MPI_COMM_WORLD, &myrank);


  assert(N%total==0);
  FILE * gpipe;
  if(!myrank)
  {
      gpipe = fopen("gnuplot.data","w");
           if(!gpipe)
           {
             exit(-1);
           }
    fprintf(gpipe,"set view 0,0\n"
	         "set pm3d at b\n"
	          "set palette rgbformulae 30,31,32\n"
               "set dgrid3d %d,%d\n"
	    "unset key\n",N>100?100:N,N>100?100:N);
  }
    // Подготовка исх. данных (только root)

  if(!myrank )
    {
    Out = (float_t *) malloc (sizeof(float_t)*N*N);
    for(i=0;i<N*N;i++)
      Out[i]=-5;
    }

  if(!myrank ||myrank==total-1){
    n = total>1 ? (int)N/total+1  :  N;
  }
  else{
    n = (int)N/total+2;
  }

  A = (float_t *) malloc(sizeof(float_t)*N*n);

  f= myrank ? myrank*N/total-1 : 0;
  printf("Total=%d, rank=%d line %d\n",total,myrank,f);
    // Инициализация матрицы A
    for (i=0; i<n; i++)
      {
	if(i+f==0)
	  for (j=0; j<N; j++)
	      A[j] = 150;

	else if(i+f==N-1)
	  for(j=0;j<N;j++)
	    A[i*N+j]=100;

	else
	  for(j=0;j<N;j++)
	  A[i*N+j]=-10;
      }


    for(i=0;i<n;i++){
      A[i*N]=-30;
      A[(i+1)*N-1]=-30;
    }


    //Plate variables
    float_t width=0.10;
    float_t height = 0.10;
    float_t lambda= 45;
    float_t C_t =460;
    float_t p=7800;
    float_t a_t=lambda/(C_t*p);

    //method variables control panel
    float_t t=0; //start/end time
    float_t dx =width/N;
    float_t dy = height/N;
    float_t dt=0.5*dy*dy/(2*a_t); //step in time
    dt=dt>0.2?0.2:dt;
    float_t dTout=2,Tout=1; // visualization result period

    float_t w_x=dt*a_t/(dx*dx);
    float_t w_y=dt*a_t/(dy*dy);
    //index of finete difference elements u-up, r-right etc.
    int u,d,l,r;

    // printf("%f %f %f %f %f %f %f %d %d\n\n",a_t,dx,dy,dt,dTout,w_x,w_y,n,N);


    while(t<Tmax)
    {
      while(t<Tout)
	{

	  for (i=1; i<n-1; i++)
	    for (j=1,k=i*N+j,l=k-1,r=k+1,u=k+N,d=k-N; j<N-1; j++,k++,l++,r++,u++,d++)
	      A[k]+=w_x*(A[r]-2*A[k]+A[l])+w_y*(A[u]-2*A[k]+A[d]);



	  //dT itaration synchronization
	  MPI_Status *status=(MPI_Status *)malloc(sizeof(MPI_Status));
	  if(myrank%2==0)
	    {
	      if(myrank!=total-1)
		MPI_Send(&(A[(n-2)*N]),N,floatMPI,myrank+1,0,MPI_COMM_WORLD);
	      if(myrank)
		MPI_Send(&(A[N]),N,floatMPI,myrank-1,0,MPI_COMM_WORLD);
	      if(myrank)
		MPI_Recv(A,N,floatMPI,myrank-1,0,MPI_COMM_WORLD,status);
	      if(myrank!=total-1)
		MPI_Recv(&(A[(n-1)*N]),N,floatMPI,myrank+1,0,MPI_COMM_WORLD,status);
	    }
	  else
	    {
	      if(myrank)
		MPI_Recv(A,N,floatMPI,myrank-1,0,MPI_COMM_WORLD,status);
	      if(myrank!=total-1)
		MPI_Recv(&(A[(n-1)*N]),N,floatMPI,myrank+1,0,MPI_COMM_WORLD,status);

	      if(myrank!=total-1)
		MPI_Send(&(A[(n-2)*N]),N,floatMPI,myrank+1,0,MPI_COMM_WORLD);
	      if(myrank)
		MPI_Send(&(A[N]),N,floatMPI,myrank-1,0,MPI_COMM_WORLD);
	    }
	  //next time
	  t+=dt;
	}
      //next output moment
      Tout+=dTout;

#ifdef GNUPLOT_ITER
      if(!myrank)
      plotresult(A,Out, N, myrank, total, gpipe);
    else
      plotresult(&(A[N]),Out,N,myrank,total,gpipe);
#endif //Output

    }
#ifdef GNUPLOT_ONE
    if(!myrank)
      plotresult(A,Out, N, myrank, total, gpipe);
    else
      plotresult(&(A[N]),Out,N,myrank,total,gpipe);
#endif
  MPI_Finalize();
  if(!myrank){
  gettimeofday(&endTime,NULL);
  float_t timewaste=(endTime.tv_sec-startTime.tv_sec)+(endTime.tv_usec-startTime.tv_usec)/1e6;
  printf("Success!T=%d\n %f sec spent\n",Tmax,timewaste);
  fclose(gpipe);
   }
  exit(0);
  }

void plotresult(float_t *A,float_t *Out, int N,int rank,int total, FILE* gnuplot)
{
int i,j;
  MPI_Gather(A,(int)N*N/total,floatMPI,
	     Out,(int)N*N/total,floatMPI,0,MPI_COMM_WORLD);
  if(!rank)
    {
      fprintf(gnuplot,"splot '-' with lines palette\n");
      for(i=0;i<N;i++)
	{
	  for(j=0;j<N;j++)
	    fprintf(gnuplot,"%d %d %f\n", i, j, Out[i * M+j]);
	}
      fprintf(gnuplot,"e\n");
      fprintf(gnuplot,"pause -1\n");
      fflush(gnuplot);
    }
  if(0){
    for(i=0;i<N;i++){
      for(j=0;j<N;j++)
	printf("%6.1f",Out[i*N+j]);
      printf("\n");
    }
  }
#ifdef DEBUG
  printf("Barrier %d\n",rank);
#endif
  MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG
  if(!rank)
    printf("Barrier passed\n");
#endif
}
