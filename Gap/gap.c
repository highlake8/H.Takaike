
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "cpu_time.c"

#define FSCANF fscanf


/***** default values of parameters ******************************************/
#define	TIMELIM	300	/* the time limit for the algorithm in seconds */
#define	GIVESOL	0	/* 1: input a solution; 0: do not give a solution */

typedef struct {
  int		timelim;	/* the time limit for the algorithm in secs. */
  int		givesol;	/* give a solution (1) or not (0) */
  /* Never modify the above two lines.  */
  /* You can add more components below. */
} Param;			/* parameters */

typedef struct {
  int	n;	/* number of jobs */
  int	m;	/* number of agents */
  int	**c;	/* cost matrix c_{ij} */
  int	**a;	/* resource requirement matrix a_{ij} */
  int	*b;	/* available amount b_i of resource for each agent i */
} GAPdata;	/* data of the generalized assignment problem */

typedef struct {
  double	timebrid;	/* the time before reading the instance data */
  double	starttime;	/* the time the search started */
  double	endtime;	/* the time the search ended */
  int		*bestsol;	/* the best solution found so far */
  /* Never modify the above four lines. */
  /* You can add more components below. */
} Vdata;		/* various data often necessary during the search */

/*************************** functions ***************************************/
void copy_parameters(int argc, char *arcv[], Param *param);
void read_instance(GAPdata *gapdata);
void prepare_memory(Vdata *vdata, GAPdata *gapdata);
void copy_instance(GAPdata *gapdata_o, GAPdata *gapdata_c);
void free_memory(Vdata *vdata, GAPdata *gapdata_o, GAPdata *gapdata_c);
void read_sol(Vdata *vdata, GAPdata *gapdata);
void recompute_cost(Vdata *vdata, GAPdata *gapdata);
void *malloc_e(size_t size);
void my_algorithm(Vdata *vdata, GAPdata *gapdata, Param *param);

/***** check the feasibility and recompute the cost **************************/
void recompute_cost(Vdata *vdata, GAPdata *gapdata)
{
  int	i, j;		
  int	*rest_b;	
  int	cost, penal;	
  int	temp;	

  rest_b = (int *) malloc_e(gapdata->m * sizeof(int));
  cost = penal = 0;
  for(i=0; i<gapdata->m; i++){rest_b[i] = gapdata->b[i];}
  for(j=0; j<gapdata->n; j++){
    rest_b[vdata->bestsol[j]] -= gapdata->a[vdata->bestsol[j]][j];
    cost += gapdata->c[vdata->bestsol[j]][j];
  }
  for(i=0; i<gapdata->m; i++){
    temp = rest_b[i];
    if(temp<0){penal -= temp;}
  }
  printf("recomputed cost = %d\n", cost);
  if(penal>0){
    printf("INFEASIBLE!!\n");
    printf(" resource left:");
    for(i=0; i<gapdata->m; i++){printf(" %3d", rest_b[i]);}
    printf("\n");
  }
  printf("time for the search:       %7.2f seconds\n",
	 vdata->endtime - vdata->starttime);
  printf("time to read the instance: %7.2f seconds\n",
	 vdata->starttime - vdata->timebrid);

  free((void *) rest_b);
}

/***** read a solution from STDIN ********************************************/
void read_sol(Vdata *vdata, GAPdata *gapdata)
{
  int	j;		
  int	value_read;	
  FILE	*fp=stdin;

  for(j=0; j<gapdata->n; j++){
    FSCANF(fp, "%d", &value_read);
    vdata->bestsol[j] = value_read - 1;
  }
}

/***** prepare memory space **************************************************/
void prepare_memory(Vdata *vdata, GAPdata *gapdata)
{
  int j;

  vdata->bestsol = (int *)  malloc_e(gapdata->n * sizeof(int));
  for(j=0; j<gapdata->n; j++){vdata->bestsol[j] = 0;}
}

/***** free memory space *****************************************************/
void free_memory(Vdata *vdata, GAPdata *gapdata_o, GAPdata *gapdata_c)
{
  free((void *) vdata->bestsol);
  free((void *) gapdata_o->c[0]);
  free((void *) gapdata_o->c);
  free((void *) gapdata_o->a[0]);
  free((void *) gapdata_o->a);
  free((void *) gapdata_o->b);
  free((void *) gapdata_c->c[0]);
  free((void *) gapdata_c->c);
  free((void *) gapdata_c->a[0]);
  free((void *) gapdata_c->a);
  free((void *) gapdata_c->b);
}


/***** copy the instance data ************************************************/
void copy_instance(GAPdata *gapdata_o, GAPdata *gapdata_c) {
  int	i, j;

  gapdata_c->m = gapdata_o->m;
  gapdata_c->n = gapdata_o->n;

  gapdata_c->c    = (int **) malloc_e(gapdata_c->m * sizeof(int *));
  gapdata_c->c[0] = (int *)  malloc_e(gapdata_c->m * gapdata_c->n * sizeof(int));
  for(i=1; i<gapdata_c->m; i++){gapdata_c->c[i] = gapdata_c->c[i-1] + gapdata_c->n;}
  gapdata_c->a    = (int **) malloc_e(gapdata_c->m * sizeof(int *));
  gapdata_c->a[0] = (int *)  malloc_e(gapdata_c->m * gapdata_c->n * sizeof(int));
  for(i=1; i<gapdata_c->m; i++){gapdata_c->a[i] = gapdata_c->a[i-1] + gapdata_c->n;}
  gapdata_c->b    = (int *)  malloc_e(gapdata_c->m * sizeof(int));

  for(i=0; i<gapdata_o->m; i++){    
    for(j=0; j<gapdata_o->n; j++){
      gapdata_c->c[i][j] = gapdata_o->c[i][j];
    }
  }

    for(i=0; i<gapdata_o->m; i++){
    for(j=0; j<gapdata_o->n; j++){
      gapdata_c->a[i][j] = gapdata_o->a[i][j];
    }
  }

  for(i=0; i<gapdata_o->m; i++){    
    gapdata_c->b[i] = gapdata_o->b[i];
  }

}

/***** read the instance data ************************************************/
void read_instance(GAPdata *gapdata)
{
  int	i, j;		
  int	value_read;	
  FILE	*fp=stdin;	

  FSCANF(fp, "%d", &value_read);	/* number of agents */
  gapdata->m = value_read;
  FSCANF(fp,"%d",&value_read);		/* number of jobs */
  gapdata->n = value_read;

  gapdata->c    = (int **) malloc_e(gapdata->m * sizeof(int *));
  gapdata->c[0] = (int *)  malloc_e(gapdata->m * gapdata->n * sizeof(int));
  for(i=1; i<gapdata->m; i++){gapdata->c[i] = gapdata->c[i-1] + gapdata->n;}
  gapdata->a    = (int **) malloc_e(gapdata->m * sizeof(int *));
  gapdata->a[0] = (int *)  malloc_e(gapdata->m * gapdata->n * sizeof(int));
  for(i=1; i<gapdata->m; i++){gapdata->a[i] = gapdata->a[i-1] + gapdata->n;}
  gapdata->b    = (int *)  malloc_e(gapdata->m * sizeof(int));
  
  for(i=0; i<gapdata->m; i++){    
    for(j=0; j<gapdata->n; j++){
      FSCANF(fp, "%d", &value_read);
      gapdata->c[i][j] = value_read;
    }
  }

  for(i=0; i<gapdata->m; i++){
    for(j=0; j<gapdata->n; j++){
      FSCANF(fp, "%d", &value_read);
      gapdata->a[i][j] = value_read;
    }
  }

  for(i=0; i<gapdata->m; i++){    
    FSCANF(fp,"%d", &value_read);
    gapdata->b[i] = value_read;
  }
}

/***** copy and read the parameters ******************************************/
void copy_parameters(int argc, char *argv[], Param *param)
{
  int i;
  param->timelim = TIMELIM;
  param->givesol = GIVESOL;

  if(argc>0 && (argc % 2)==0){
    printf("USAGE: ./gap [param_name, param_value] [name, value]...\n");
    exit(EXIT_FAILURE);}
  else{
    for(i=1; i<argc; i+=2){
      if(strcmp(argv[i],"timelim")==0) param->timelim = atoi(argv[i+1]);
      if(strcmp(argv[i],"givesol")==0) param->givesol = atoi(argv[i+1]);
    }
  }
}

void *malloc_e( size_t size ) {
  void *s;
  if ( (s=malloc(size)) == NULL ) {
    fprintf( stderr, "malloc : Not enough memory.\n" );
    exit( EXIT_FAILURE );
  }
  return s;
}

/***** main Algorithm *********************************************/
void swap(int i, int j, double *A)
{
 double temp;
 temp=A[i]; A[i]=A[j]; A[j]=temp;
 return;
}

void swapint(int i, int j, int *A)
{
 int temp;
 temp=A[i]; A[i]=A[j]; A[j]=temp;
 return;
}

int pivot(int i, int j, double *A)     
{
    int pv, k;
    
    k=i+1;
    while(k<=j && A[i]==A[k]) k=k+1;
    if(k>j) pv=-1;
    else if(A[i]>=A[k]) pv=i;
    else pv=k;
    
    return(pv);
}

int partition(int i, int j, double a, double *A)
{
 int l, r, k;

 l=i; r=j;                     
 while(1)    
   {
    while(A[l]<a){
      l=l+1;
    } 
    while(A[r]>=a){
      r=r-1;
    }  
    if(l<=r){
      swap(l, r, A); l=l+1; r=r-1;
    }    
    else break;                                
   }
 k=l;
 return(k);
}

void quicksort(int i, int j, double *A)
{
 int pv, k;
 double a;

 pv=pivot(i, j, A);
 if(pv!=-1)                    
   {
    a=A[pv];                   
    k=partition(i, j, a, A);   
    quicksort(i, k-1, A);      
    quicksort(k, j, A);        
   }
 return;
}

int random3(int t)
{
  int a;
  t = (21*t+7) % 1000;
  a = t % 3;

  return a;
}

void randamized_greedy(int *x, int **a, int *b, int **c, int m, int n, int r)
{
  double **d, **e;

  d = (double **)calloc(n, sizeof(double *));
  for(int i = 0; i < n; i++) {
	  d[i] = (double *)calloc(m, sizeof(double));
  }
  e = (double **)calloc(n, sizeof(double *));
  for(int i = 0; i < n; i++) {
	  e[i] = (double *)calloc(m, sizeof(double));
  }
  for(int i = 0; i < n; i++){
    x[i] = -1;
  }

  int i, j, k, l, s;

  j = 0;
  while(j < n){
    for(i = 0; i < m; i++){  
      d[j][i] = -(b[i]) + (a[i][j]);
      e[j][i] = d[j][i];
    }
    quicksort(0, m-1, e[j]);

    l = 0;
    while(l < m && e[j][l] >= 0){
      l++;
    }

    s = random3(r);
    r = s;
    l = l+s;
    int v = 0;

    if(l + s >= m){
      v = -m + l + s;
    }

    k = 0;
    while(d[j][k] != e[j][l-v]){
      k++;
    }
    
    b[k] = b[k] - a[k][j];
    x[j] = k;
    
    
    int p = 0;
    int t = 0;
    int o, u;
    int w = 0;
    
    while(w < m){
      if(b[w] >= 0){w=w+1; continue;}
      else{
        t = w;
        o = 0;
        while(t < m && b[t] < 0){
          while(o < n && x[o] != t){
            o++;
            if(o >= n){o = 0; t = t + 1; continue;}
          }



          u = -m + 1;
          while(t+u < 0 || e[o][t+u] >= 0){
            u = u + 1;
            if(t + u >= m)break;
          }
          if(t + u >= m)continue;
          p = 0;
          while(p < m && d[o][p] != e[o][t+u]){
            p = p + 1;
          }       
          b[t] = b[t] + a[t][o];   
          x[o] = p;  
          b[p] = b[p] - a[p][o];    
          if(b[p] < 0){
            b[p] = b[p] + a[p][o];
            b[t] = b[t] - a[t][o];
            x[o] = t; 
            o = o + 1;
            continue;
          }
          else t = t + 1;
        }
  
        if(b[t] < 0){w = w + 1; continue;}
        
      }
    }

    j = j + 1;
  }
  
  for(int j=0;j<n;j++){
    free(d[j]);
  }
  free(d);
  for(int j=0;j<n;j++){
    free(e[j]);
  }
  free(e);
}

void evaluation(int costsum, int costsum2, int n, int **a, int *b, int **c, int *x)
{
  for(int t = 0; t < n; t++){
    if(x[t] == -1){
      return;
    }
  }
  do{
    costsum = 0;
    for(int l = 0; l < n; l++){
      costsum += c[x[l]][l];
    }
    for(int i = 0; i < n; i++){
      for(int j = 0; j < n; j++){
        costsum2 = costsum - c[x[i]][i] - c[x[j]][j] + c[x[i]][j] + c[x[j]][i];

        if(costsum2 < costsum && b[x[i]] + a[x[i]][i] - a[x[i]][j] >= 0 && b[x[j]] + a[x[j]][j] - a[x[j]][i] >= 0){
          b[x[i]] = b[x[i]] + a[x[i]][i] - a[x[i]][j];
          b[x[j]] = b[x[j]] + a[x[j]][j] - a[x[j]][i];
          swapint(i, j, x);
          break;
        }
        else {
          costsum2 = costsum;
        }
      }
      if(costsum2 < costsum){
        break;
      }
    }
  }while(costsum2 < costsum);
  if(costsum2 >= costsum){
    return;
  }
}

void simple_local_search(int *x, int **a, int *b, int **c, int n)
{
  int costsum = 0;
  int costsum2 = 0;
  
  evaluation(costsum, costsum2, n, a, b, c, x);

}

void my_algorithm(Vdata *vdata, GAPdata *gapdata, Param *param){

  int t = 0;
  int *x, *x2, *b;
  int costsum = 0, costsum2 = 0;

  x = malloc(gapdata->n*sizeof(int));
  x2 = malloc(gapdata->n*sizeof(int));
  b = malloc(gapdata->m*sizeof(int));


  for(int i = 0; i < gapdata->n; i++){
    x[i] = 0;
  }
  for(int i = 0; i < gapdata->n; i++){
    x2[i] = 0;
  }
  for(int i = 0; i < gapdata->m; i++){
    b[i] = gapdata->b[i];
  }

  for(int i = 0; i < gapdata->n; i++){
    vdata->bestsol[i] = 1;
  }

  randamized_greedy(x, gapdata->a, b, gapdata->c, gapdata->m, gapdata->n, t);
  simple_local_search(x, gapdata->a, b, gapdata->c, gapdata->n);
  for(int t = 0; t < gapdata->n; t++){
    if(x[t] == -1){
      costsum = 1000000000;
    }
  }

  for(int i = 0; i < gapdata->n; i++){
    costsum += gapdata->c[x[i]][i];
  }

  while((cpu_time() - vdata->starttime) < param->timelim){
    for(int i = 0; i < gapdata->m; i++){
      b[i] = gapdata->b[i];
    }
    randamized_greedy(x2, gapdata->a, b, gapdata->c, gapdata->m, gapdata->n, t);
    simple_local_search(x2, gapdata->a, b, gapdata->c, gapdata->n);
    for(int t = 0; t < gapdata->n; t++){
      if(x[t] == -1){
        costsum2 = 1000000000;
      }
      else{
        costsum2 = 0;
        costsum2 += gapdata->c[x2[t]][t];
      }
    } 
    if(costsum > costsum2){
      for(int k = 0; k < gapdata->n; k++){
        vdata->bestsol[k] = x2[k];
      }
      costsum = costsum2;
      for(int k = 0; k < gapdata->n; k++){
        x[k] = x2[k];
      }
    }
    else{
      for(int k = 0; k < gapdata->n; k++){
        vdata->bestsol[k] = x[k];
      }
    }
  }

  free(x);
  free(x2);
  free(b);

  return;
}
  

int main(int argc, char *argv[])
{
  cpu_time();
  Param		param;		
  GAPdata	gapdata_org;	
  GAPdata	gapdata_cpd;	
  Vdata		vdata;	

  vdata.timebrid = cpu_time();
  copy_parameters(argc, argv, &param);
  read_instance(&gapdata_org);
  copy_instance(&gapdata_org, &gapdata_cpd);
  prepare_memory(&vdata, &gapdata_cpd);
  if(param.givesol==1){read_sol(&vdata, &gapdata_cpd);}
  vdata.starttime = cpu_time();
  my_algorithm(&vdata, &gapdata_cpd, &param);
  vdata.endtime = cpu_time();
  recompute_cost(&vdata, &gapdata_org);
  free_memory(&vdata, &gapdata_org, &gapdata_cpd);

  return EXIT_SUCCESS;
}
