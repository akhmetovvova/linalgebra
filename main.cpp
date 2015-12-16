/*
 * main.cpp
 *
 *  Created on: Dec 16, 2015
 *      Author: vova
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <limits.h>
#include <errno.h>
#include <math.h>
#include <unistd.h>
#include <gmp.h>
#include "gmpxx.h"
#include "ccarray.h"
#include "Matrix.h"
//#include "Vector.h"
//#include "obj.h"
//#include "linalg.h"



/** object data */
typedef
struct obj_t {
  double x, y, f, ef,mag;
} obj_t;


typedef _Matrix <mpf_class> Matrix;

#define UNUSED(x)               ((void)(x))
#define MAX_HEADER_LENGTH       2048
#define MAX_INPUT_LINE_LENGTH   2048
#define PI                      M_PI

pthread_mutex_t lock; //Исключающая блокировка


typedef
enum compression_t {
    compression_unknown = -1,
    compression_none,
    compression_bzip,
    compression_gzip
} compression_t;
//using namespace std;

/** uses fopen() if input file seems to be uncompresssed, and popen() if input file seems compressed */
static FILE * open_file( const char * fname, compression_t * compression )
{
  char cmd[PATH_MAX] = {0};
  FILE * input = NULL;

  if ( *compression == compression_unknown )
  {
    const char * suffix;

    if ( (suffix = strstr(fname, ".bz2")) && *(suffix + 4) == 0 ) {
      *compression = compression_bzip;
    }
    else if ( (suffix = strstr(fname, ".bz")) && *(suffix + 3) == 0 ) {
      *compression = compression_bzip;
    }
    else if ( (suffix = strstr(fname, ".gz")) && *(suffix + 3) == 0 ) {
      *compression = compression_gzip;
    }
    else {
      *compression = compression_none;
    }
  }

  switch ( *compression )
  {
  case compression_bzip:
    snprintf(cmd, sizeof(cmd) - 1, "bzip2 -dc '%s'", fname);
    if ( !(input = popen(cmd, "r"))) {
      fprintf(stderr, "popen('%s') fails: %s\n", cmd, strerror(errno));
    }
    break;

  case compression_gzip:
    snprintf(cmd, sizeof(cmd) - 1, "gzip -dc '%s'", fname);
    if ( !(input = popen(cmd, "r"))) {
      fprintf(stderr, "popen('%s') fails: %s\n", cmd, strerror(errno));
    }
    break;

  case compression_none:
    if ( !(input = fopen(fname, "r")) ) {
      fprintf(stderr, "fopen('%s') fails: %s\n", fname, strerror(errno));
    }
    break;

  default:
    fprintf(stderr, "BUG IN CODE: invalid compression tag=%d\n", *compression);
    break;
  }

  return input;
}

/** uses fclose() if input file was uncompressed, and pclose() for compressed case */
static void close_file( FILE * input, compression_t compression )
{
  if ( input && input != stdin ) {
    if ( compression > compression_none ) {
      pclose(input);
    }
    else {
      fclose(input);
    }
  }
}

void skipline(FILE * fr)
{
	while(!feof(fr))
		if(fgetc(fr)=='\n')
			break;
}

mpf_class pow(const mpf_class& op1, unsigned long int op2, int bitsize)
{
	mpf_set_default_prec(bitsize);

	mpf_class ret;
	mpf_pow_ui(ret.__get_mp(),op1.__get_mp(),op2);
	return ret;
}

static long int correction_objects( FILE * output, int x, int y, int f, int k, int p, mpf_class * A,int bitsize, int print_df)
{
	mpf_set_default_prec(bitsize);

	char line[MAX_INPUT_LINE_LENGTH];
	int ic;
	char * pc;
	double xx,yy,zz;
	mpf_class mx=0;
	mpf_class my=0;
	mpf_class ff=0;
	mpf_class xy[k];

	char header[MAX_HEADER_LENGTH];

	  if ( !fgets(header, MAX_HEADER_LENGTH, output) ) {
	    fprintf(stderr,"fgets(header) fails: %d (%s)\n", errno, strerror(errno));
	    return -1;
	  }

	  if ( header[MAX_HEADER_LENGTH - 1] != 0 ) {
	    fprintf(stderr,"too long header line in this file\n");
	    return -1;
	  }


	int ii=0;

	if(print_df)
	{
	printf("cf\tdf\t%s",header);
	}else{printf("cf\t%s",header);}


	while ( !feof(output) )
	{

	  //printf("%.10f\t%.10f\t%.2f\t%.2f\n",m[i].x,m[i].y,m[i].f,m[i].ef);

      line[MAX_INPUT_LINE_LENGTH-1] = 0;
      if ( !fgets(line, MAX_INPUT_LINE_LENGTH, output) ) {
        break;
      }


      if ( line[MAX_INPUT_LINE_LENGTH - 1] != 0 ) {
        fprintf(stderr,"too long input line in this file\n");
        return -1;
      }

      /* remove trailing new line */
      line[strlen(line) - 1] = 0;


      for ( pc = line, ic = 1; ic < x; ++ic ) {
        if ( !(pc = strchr(pc + 1, '\t')) ) {
          break;
        }
      }
      if ( ic != x || sscanf(pc, " %lf",&xx) != 1 ) {
        continue;
      }

      for ( pc = line, ic = 1; ic < y; ++ic ) {
        if ( !(pc = strchr(pc + 1, '\t')) ) {
          break;
        }
      }
      if ( ic != y || sscanf(pc, " %lf", &yy) != 1 ) {
        continue;
      }
      for ( pc = line, ic = 1; ic < f; ++ic ) {
		  if ( !(pc = strchr(pc + 1, '\t')) ) {
			break;
		  }
		}
		if ( ic != f || sscanf(pc, " %lf", &zz) != 1 ) {
		  continue;
		}


      int l=0;
      for(int z=0;z<=p;++z)
      {
      	for(int n=z,m=0;m<=z;n--,m++)
      		{
      			mx=xx;
      			mx=pow(mx,n,bitsize);

      			my=yy;
      			my=pow(my,m,bitsize);

      			xy[l]=mx*my;

      			 //gmp_printf ("x= %Ff \n",mx.__get_mp());
      			 //gmp_printf ("y= %Ff \n",my.__get_mp());
      			 //gmp_printf ("x[%d]= %Ff \n",l,f[l].__get_mp());

      			l++;
      		}
      }


      ff=0;
      for(int i=0;i<k;++i)
      {
    	  ff+=A[i]*xy[i];
      }

      //gmp_printf ("ff=%Ff\n",ff.__get_mp());
      if(print_df)
      {
    	  mpf_class df=ff-zz;
    	  gmp_printf("%Ff\t%Ff\t%s\n",ff.__get_mp(),df.__get_mp(),strdup(line));
      }else{
      	  gmp_printf("%Ff\t%s\n",ff.__get_mp(),strdup(line));
      }

      ii++;
	}


  return ii;
}
static long int correction_objects( FILE * output, int x, int y, int f, int f1, int k, int p, mpf_class * A, mpf_class * A1,int bitsize,int print_df)
{
	mpf_set_default_prec(bitsize);

	char line[MAX_INPUT_LINE_LENGTH];
	int ic;
	char * pc;
	double xx,yy,zz,zz1;
	mpf_class mx=0;
	mpf_class my=0;
	mpf_class ffx,ffy;
	mpf_class xy[k];

	char header[MAX_HEADER_LENGTH];

	  if ( !fgets(header, MAX_HEADER_LENGTH, output) ) {
	    fprintf(stderr,"fgets(header) fails: %d (%s)\n", errno, strerror(errno));
	    return -1;
	  }

	  if ( header[MAX_HEADER_LENGTH - 1] != 0 ) {
	    fprintf(stderr,"too long header line in this file\n");
	    return -1;
	  }


	int ii=0;

	if(print_df)
	{
	printf("cx\tcy\tdf\tdf1\t%s",header);
	}else{printf("cx\tcy\t%s",header);}


	while ( !feof(output) )
	{

	  //printf("%.10f\t%.10f\t%.2f\t%.2f\n",m[i].x,m[i].y,m[i].f,m[i].ef);

      line[MAX_INPUT_LINE_LENGTH-1] = 0;
      if ( !fgets(line, MAX_INPUT_LINE_LENGTH, output) ) {
        break;
      }


      if ( line[MAX_INPUT_LINE_LENGTH - 1] != 0 ) {
        fprintf(stderr,"too long input line in this file\n");
        return -1;
      }

      /* remove trailing new line */
      line[strlen(line) - 1] = 0;


      for ( pc = line, ic = 1; ic < x; ++ic ) {
        if ( !(pc = strchr(pc + 1, '\t')) ) {
          break;
        }
      }
      if ( ic != x || sscanf(pc, " %lf",&xx) != 1 ) {
        continue;
      }

      for ( pc = line, ic = 1; ic < y; ++ic ) {
        if ( !(pc = strchr(pc + 1, '\t')) ) {
          break;
        }
      }
      if ( ic != y || sscanf(pc, " %lf", &yy) != 1 ) {
        continue;
      }
      for ( pc = line, ic = 1; ic < f; ++ic ) {
		  if ( !(pc = strchr(pc + 1, '\t')) ) {
			break;
		  }
		}
		if ( ic != f || sscanf(pc, " %lf", &zz) != 1 ) {
		  continue;
			}
		  for ( pc = line, ic = 1; ic < f1; ++ic ) {
		  if ( !(pc = strchr(pc + 1, '\t')) ) {
			break;
		  }
		}
		if ( ic != f1 || sscanf(pc, " %lf", &zz1) != 1 ) {
		  continue;
		}

      int l=0;
      for(int z=0;z<=p;++z)
      {
      	for(int n=z,m=0;m<=z;n--,m++)
      		{
      			mx=xx;
      			mx=pow(mx,n,bitsize);

      			my=yy;
      			my=pow(my,m,bitsize);

      			xy[l]=mx*my;

      			 //gmp_printf ("x= %Ff \n",mx.__get_mp());
      			 //gmp_printf ("y= %Ff \n",my.__get_mp());
      			 //gmp_printf ("x[%d]= %Ff \n",l,f[l].__get_mp());

      			l++;
      		}
      }


      ffx=0;
      ffy=0;
      for(int i=0;i<k;++i)
      {
    	  ffx+=A[i]*xy[i];
    	  ffy+=A1[i]*xy[i];
      }
      if(print_df)
      {
    	  mpf_class df =zz-ffx;
    	  mpf_class df1=zz1-ffy;
      	  gmp_printf("%Ff\t%Ff\t%Ff\t%Ff\t%s\n",ffx.__get_mp(),ffy.__get_mp(),df.__get_mp(),df1.__get_mp(),strdup(line));
      }else{gmp_printf("%Ff\t%Ff\t%s\n",ffx.__get_mp(),ffy.__get_mp(),strdup(line));}

      ii++;
	}

return ii;
}

static long int load_objects_xy(FILE * input, int x, int y, int f, int ef, double *avg,obj_t * m)
{
	double avg_ef=0;
	char line[MAX_INPUT_LINE_LENGTH];
	int ic;
	char * pc;
	int i=0;

	while ( !feof(input) )
	{
      line[MAX_INPUT_LINE_LENGTH-1] = 0;
      if ( !fgets(line, MAX_INPUT_LINE_LENGTH, input) ) {
        break;
      }

      if ( line[MAX_INPUT_LINE_LENGTH - 1] != 0 ) {
        fprintf(stderr,"too long input line in this file\n");
        return -1;
      }

      line[strlen(line) - 1] = 0;

      for ( pc = line, ic = 1; ic < x; ++ic ) {
        if ( !(pc = strchr(pc + 1, '\t')) ) {
          break;
        }
      }
      if ( ic != x || sscanf(pc, " %lf",&m[i].x) != 1 ) {
        continue;
      }

      for ( pc = line, ic = 1; ic < y; ++ic ) {
        if ( !(pc = strchr(pc + 1, '\t')) ) {
          break;
        }
      }
      if(y>0)
      {
      if ( ic != y || sscanf(pc, " %lf", &m[i].y) != 1 ) {
        continue;
      }

      for ( pc = line, ic = 1; ic < f; ++ic ) {
  	  if ( !(pc = strchr(pc + 1, '\t')) ) {
  		break;
  	  	  }
      }
      }
  	if ( ic != f || sscanf(pc, " %lf", &m[i].f) != 1 ) {
  	  continue;
  	}

  	if(ef>0)
  	{
  	for ( pc = line, ic = 1; ic < ef; ++ic ) {
  	  if ( !(pc = strchr(pc + 1, '\t')) ) {
  		break;
  	  }
  	}
  	if ( ic != ef || sscanf(pc, " %lf", &m[i].ef) != 1 ) {
  	  continue;
  	}
  	avg_ef+=(m[i].ef-avg_ef)/(i+1);
  	}

  	i++;

	}

	*avg=avg_ef;
  return i;
}
static long int load_objects_xym(FILE * input, int x, int y, int mag, int f, int ef, double *avg, obj_t * m)
{
	double avg_ef=0;
	char line[MAX_INPUT_LINE_LENGTH];
	int ic;
	char * pc;
	int i=0;

	while ( !feof(input) )
	{
      line[MAX_INPUT_LINE_LENGTH-1] = 0;
      if ( !fgets(line, MAX_INPUT_LINE_LENGTH, input) ) {
        break;
      }

      if ( line[MAX_INPUT_LINE_LENGTH - 1] != 0 ) {
        fprintf(stderr,"too long input line in this file\n");
        return -1;
      }

      line[strlen(line) - 1] = 0;

      for ( pc = line, ic = 1; ic < x; ++ic ) {
        if ( !(pc = strchr(pc + 1, '\t')) ) {
          break;
        }
      }
      if ( ic != x || sscanf(pc, " %lf",&m[i].x) != 1 ) {
        continue;
      }

      if(y>0)
      {

          for ( pc = line, ic = 1; ic < y; ++ic ) {
            if ( !(pc = strchr(pc + 1, '\t')) ) {
              break;
            }
          }

		  if ( ic != y || sscanf(pc, " %lf", &m[i].y) != 1 ) {
			continue;
		  }

      }
      if(mag>0)
      {
      for ( pc = line, ic = 1; ic < mag; ++ic ) {
        if ( !(pc = strchr(pc + 1, '\t')) ) {
          break;
        }
      }

      if ( ic != mag || sscanf(pc, " %lf", &m[i].mag) != 1 ) {
        continue;
      }
      }

      for ( pc = line, ic = 1; ic < f; ++ic ) {
  	  if ( !(pc = strchr(pc + 1, '\t')) ) {
  		break;
  	  	  }
      }

  	if ( ic != f || sscanf(pc, " %lf", &m[i].f) != 1 ) {
  	  continue;
  	}

  	if(ef>0)
  	{
  	for ( pc = line, ic = 1; ic < ef; ++ic ) {
  	  if ( !(pc = strchr(pc + 1, '\t')) ) {
  		break;
  	  }
  	}
  	if ( ic != ef || sscanf(pc, " %lf", &m[i].ef) != 1 ) {
  	  continue;
  	}
  	avg_ef+=(m[i].ef-avg_ef)/(i+1);
  	}

  	i++;

	}

	*avg=avg_ef;
  return i;
}

struct int_int
{int i;
 int j;
};

int get_power(int i,int p, int *nn,int *mm)
{
	if(i<0)
	{return -1;
	}
	else{

		int l=0;
		for(int z=0;z<=p;++z)
			for(int n=z,m=0;m<=z;n--,m++)
			{
				if(i==l)
				{
					*nn=n;
					*mm=m;
				return 1;
				}

				l++;
			}
	}
	return -1;
}


void  calc_sigma(mpf_class *sigma, double * wn, mpf_class * L, obj_t * mas, int k, long int nstars, double avg_ef, int p, int NBIT, int * N)
{
	//printf("Calc sigma k=%d\n",k);

	//gmp_printf("carent sigma = %Ff\n",sigma->__get_mp());

	double weightN=0;

	mpf_set_default_prec(NBIT);

	mpf_class xx,yy,ww,ll,ff,ws,sigm,zz;
	ff=0;
	ll=0;
	zz=0;
	ww=0;
	xx=0;
	yy=0;
	ws=0;
	sigm=0;
	double weight=0;

	mpf_class x[k];

	for(int i=0;i<k;++i)
	{
	x[i]=0;
	}

	//int s=0;
	for(int s=0;s<nstars;++s)
	{
			if(avg_ef!=0)
			{
				weight=avg_ef/mas[s].ef; //weight of stars
			}else{weight=1;}

			ww=weight;
			//gmp_printf("weight=%.5Ff\n",ww);

			ff=0;
			for(int i=0;i<k;++i)
			{
				int n=0,m=0;
				if(get_power(i,p,&n,&m)>0)
				{
					xx=mas[s].x;
					yy=mas[s].y;

					x[i]=pow(xx,n,NBIT)*pow(yy,m,NBIT)*N[i];
					ff+=x[i]*L[i];
				}

			}
			//gmp_fprintf(stderr,"ff=%Fe\n",ff.__get_mp());

			//printf("\n");
			ll=(mas[s].f-ff);
			if(abs(ll)<(*sigma))
			{
				//gmp_fprintf(fr,"%lf\t%lf\t%lf\tlf\t%Ff\t%Ff\t%Ff\n",mas[s].x,mas[s].y,mas[s].f,mas[s].ef,ff.__get_mp(),ll.__get_mp(),ff.__get_mp()-ll.__get_mp());

				ws=pow(ll,2,NBIT);
				sigm+=ws*ww;
				//sigm+=ll*ww;
				weightN+=weight;
			}
	}

	*wn=weightN;
	*sigma=sigm;
	//gmp_printf("sigmma=%.Ff and weight=%lf\n",sigm.__get_mp(),*wn);
	//printf("weight=%lf\n",*wn);
	return;
}


long int crerat_matrix_sig(int p,int k, long int nstars, double avg_ef, obj_t * mas, mpf_class **CM, mpf_class *YM, mpf_class * A, int bitsize, mpf_class nsig, double *avgmag, int *N)
{
	//printf("k=%d\n",k);

	long int new_nstars=0;

	mpf_set_default_prec(bitsize);
	mpf_class ww=0;
	mpf_class ww_sum=0;
	mpf_class xx=0;
	mpf_class yy=0;
	mpf_class xij=0;
	mpf_class ff=0;
	mpf_class ll=0;

	double weight;

	mpf_class x[k];

	for(int i=0;i<k;++i)
	{
	x[i]=0;
	}

	//avgmag=0;

	//int s=0;
	for(int s=0;s<nstars;++s)
	{
			if(avg_ef!=0)
			{
				weight=avg_ef/mas[s].ef; //weight of stars
			}else{weight=1;}

			ww=weight;
			ww_sum+=ww;

			*avgmag+=mas[s].mag*weight;

			//gmp_printf("weight=%.5Ff\n",ww);

			int n=0,m=0;
			ff=0;
			for(int i=0;i<k;++i)
			{

				//
				//if(get_power(N[i],p,&n,&m)>0)
				if(get_power(i,p,&n,&m)>0 )
					{

					xx=mas[s].x;
					yy=mas[s].y;

					x[i]=pow(xx,n,bitsize)*pow(yy,m,bitsize)*N[i];
					ff+=x[i]*A[i];
					//gmp_printf("x[%d]=%.5Ff\n",l,x[l].__get_mp());
					}
			}

			//printf("\n");
			ll=(mas[s].f-ff);
			if(abs(ll)<nsig)
			{
				for(int i=0;i<k;i++)
				{
					YM[i]+=(mas[s].f*x[i]*ww);

						for(int j=0;j<k;j++)
						{
							xij=x[i]*x[j]*ww;
							CM[i][j]+=xij;
							//gmp_fprintf(stderr,"C=%Fe\n",CM[i][j].__get_mp());
						}
				}

				++new_nstars;
			}
	}

	//fprintf(stderr,"magbin=%lf\n",avgmag/new_nstars);
	*avgmag/=(mpf_get_d(ww_sum.__get_mp()));
	return new_nstars;
}

static int partition (double * m, int a, int b)
	{
	  int i = a;
	  int j;
	  for (j = a; j <= b; j++ )    // просматриваем с a по b
	   {
		 if (m[j] <= m[b])            // если элемент m[j] не превосходит m[b],
		  {
			//swap(m[i],m[j]);         // меняем местами m[j] и m[a], m[a+1], m[a+2] и так далее...
			double temp=m[i];						  // то есть переносим элементы меньшие m[b] в начало,
			m[i]=m[j];						  // а затем и сам m[b] «сверху»
			m[j]=temp;
			i++;                      // таким образом последний обмен: m[b] и m[i], после чего i++
		  }
	   }
	  return i-1;                     // в индексе i хранится <новая позиция элемента m[b]> + 1
	}

static void quicksort (double * m, int a, int b) // a - начало подмножества, b - конец
	{                                        // для первого вызова: a = 0, b = <элементов в массиве> - 1
	 if (a >= b) return;
	 int c = partition (m, a, b);
	 //printf("a=%d\tb=%d\tc=%d\n",a,b,c);

	 quicksort (m, a, c-1);
	 quicksort (m, c+1, b);

	}

void multiplication(mpf_class **A, mpf_class **B, mpf_class **C,int N)
{
	for (int i=0; i < N; i++)
	    for (int k=0; k < N; k++)
	        for (int j=0;j < N; j++)
	            C[i][k]+=(A[i][j])*(B[j][k]);
}

void multiplication(mpf_class *C, mpf_class **A, mpf_class *B, int N)
{
	for (int i=0; i < N; i++)
	        for (int j=0;j < N; j++)
	            C[i]+=(A[i][j])*(B[j]);
}


void inversion(mpf_class **A, int N)
{
	mpf_class temp;

	mpf_class **E = new mpf_class *[N];

    for (int i = 0; i < N; i++)
        E[i] = new mpf_class [N];

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            E[i][j] = 0.0;

            if (i == j)
                E[i][j] = 1.0;
        }

    for (int k = 0; k < N; k++)
    {
        temp = A[k][k];

        for (int j = 0; j < N; j++)
        {
            A[k][j] /= temp;
            E[k][j] /= temp;
        }

        for (int i = k + 1; i < N; i++)
        {
            temp = A[i][k];

            for (int j = 0; j < N; j++)
            {
                A[i][j] -= A[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }

    for (int k = N - 1; k > 0; k--)
    {
        for (int i = k - 1; i >= 0; i--)
        {
            temp = A[i][k];

            for (int j = 0; j < N; j++)
            {
                A[i][j] -= A[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A[i][j] = E[i][j];

    for (int i = 0; i < N; i++)
        delete [] E[i];

    delete [] E;
}

void pprint_residual(mpf_class sigma,int nsigm, double avg_ef,int k,int p,int bitsize,long int nstars,obj_t * massiv,mpf_class * L, FILE * fr)
{
		mpf_set_default_prec(bitsize);

		mpf_class ff=0;
		mpf_class ll=0;

		mpf_class xx=0;
		mpf_class yy=0;
		mpf_class ws=0;

		double weight=0;

		mpf_class x[k];

		fprintf(fr,"x\ty\tf\tef\tweight\tcalcf\tO_C\n");

		for(int i=0;i<k;++i)
		{
		x[i]=0;
		}

		//int s=0;
		for(int s=0;s<nstars;++s)
		{
				if(avg_ef!=0)
				{
					weight=avg_ef/massiv[s].ef; //weight of stars
				}else{weight=1;}

				//gmp_printf("weight=%.5Ff\n",ww);

				int l=0;
				for(int z=0;z<=p;++z)
				{
					for(int n=z,m=0;m<=z;n--,m++)
					{
						xx=massiv[s].x;
						xx=pow(xx,n,bitsize);

						yy=massiv[s].y;
						yy=pow(yy,m,bitsize);

						x[l]=xx*yy;

						 //gmp_printf ("x= %Ff \n",xx.__get_mp());
						 //gmp_printf ("y= %Ff \n",yy.__get_mp());
						 //gmp_printf ("x[%d]= %Ff \n",l,x[l].__get_mp());

						l++;
					}
				}

				ff=0;
				for(int i=0;i<k;i++)
				{
					ff+=x[i]*L[i];
				}
				//printf("\n");
				ll=(massiv[s].f-ff);
				if(abs(ll)<(sigma))
				{
					gmp_fprintf(fr,"%lf\t%lf\t%lf\t%lf\t%lf\t%Ff\t%Ff\n",massiv[s].x,massiv[s].y,massiv[s].f,massiv[s].ef,weight,ff.__get_mp(),ll.__get_mp());
				}
		}
	return;
}

//(&sigma,avg_ef,k,ci,bitsize,massiv);
int calc_polifit (mpf_class * sigma,int nsigm,double avg_ef,int k,int p,int ci,int bitsize, long int Nstars, obj_t * massiv,mpf_class * A,mpf_class * EA, int beverbose, double *magbin, int *nstarss, double student_test, int model)
{
  mpf_set_default_prec(bitsize);
  mpf_class **C 	= new mpf_class *[k];
  mpf_class **invC 	= new mpf_class *[k];
  mpf_class *Y 		= new mpf_class [k];

  int * N = new int [k];

  for (int i = 0; i < k; i++)
  {
  C[i] = new mpf_class [k];
  invC[i] = new mpf_class [k];
  A[i]=0;
  EA[i]=0;
  N[i]=1;
  }

  	mpf_class sigm=*sigma;
	mpf_class temp_sigma=2*sigm;
	double temp_weight;
  	int ici=0;
  	long int nstars=Nstars;
 	double weightN=Nstars;

  	double mag=0;

	  	while(ici<ci && sigm<temp_sigma )
	  	{

	  		//fprintf(stderr,"Itaration %d and magbin=%lf\n",ici,mag);

	  		temp_sigma=sigm;
	  		temp_weight=weightN;

	  		sigm=sqrt(sigm/(temp_weight-k));

	  		//gmp_printf("sigma=%Ff\n",sigma.__get_mp());

			for(int i=0;i<k;++i)
			{
				Y[i]=0;
				for(int j=0;j<k;++j)
				{
					C[i][j]=0;
				}
			}

			sigm*=nsigm;
			nstars=crerat_matrix_sig(p,k,Nstars,avg_ef,massiv,C,Y,A,bitsize,sigm,&mag,N);

			*nstarss=nstars;
			//printf("nstars=%ld\n",nstars);

			for(int i=0;i<k;++i)
			{
				for(int j=0;j<k;++j)
				   {
					invC[i][j]=C[i][j];
				   }
				A[i]=0;
			}

			inversion(invC,k);
			//multiplication(C,invC,E,k);


			multiplication(A,invC,Y,k);


			calc_sigma(&sigm,&weightN,A,massiv,k,Nstars,avg_ef,p,bitsize,N);


			ici++;

			double sig=sqrt(mpf_get_d(sigm.__get_mp())/(weightN-k));
			  	if(beverbose)
			  	{
				fprintf(stderr,"Nstars=%ld weight=%lf sigma= %lf k=%d avg_mag=%6.3f\n",Nstars,weightN,sig,ici,mag);
			  	}

	  	}

	  	int kk=k;
	  	if(student_test>=1 || model>0)
	  	{

	  		//point for go in while
	  		temp_sigma=sigm*3;

	  		while(ici<ci && sigm!=temp_sigma )
	  		{

	  		if(beverbose)
	  		{
	  		fprintf(stderr,"st=%lf\n",student_test);
	  		}
	  		temp_sigma=sigm;
	  		temp_weight=weightN;

	  		sigm=sqrt(sigm/(temp_weight-kk));

	  		//gmp_fprintf(stderr,"sigma=%Ff\n",sigm.__get_mp());

	  		int l=0;

	  		for(int i=0;i<k;i++)
			{
	  			if(model>0)
	  			{
	  				if(model==1)
	  				{
	  					l=9;
	  					switch (i)
	  					{
	  					case 0: N[i]=1;
	  					break;
	  					case 1: N[i]=1;
	  					break;
	  					case 2: N[i]=1;
	  					break;
	  					case 3: N[i]=1;
	  					break;
	  					case 4: N[i]=1;
	  					break;
	  					case 6: N[i]=1;
	  					break;
	  					case 8: N[i]=1;
	  					break;
	  					case 17: N[i]=1;
	  					break;
	  					case 19: N[i]=1;
	  					break;

	  					default : N[i]=0;
	  					break;
	  					}

	  				}else
	  				{
	  					l=9;
						switch (i)
						{
						case 0: N[i]=1;
						break;
						case 1: N[i]=1;
						break;
						case 2: N[i]=1;
						break;
						case 4: N[i]=1;
						break;
						case 5: N[i]=1;
						break;
						case 7: N[i]=1;
						break;
						case 9: N[i]=1;
						break;
						case 16: N[i]=1;
						break;
						case 18: N[i]=1;
						break;

						default : N[i]=0;
						break;
						}
	  				}

	  			}else{

					EA[i]=sigm*sqrt(invC[i][i]);

					if(beverbose)
					{
						gmp_fprintf(stderr,"!!!!!!EA[%d]=%Fe\n",i,EA[i].__get_mp());
					}

				  if(EA[i]>0)
				  {
					  if(abs(A[i]/EA[i]) <= student_test)
					  {
						  if(beverbose)
						  {
						  gmp_fprintf(stderr,"!!!A[%d]=%Fe\t EA[%d]=%Fe\n",i,A[i].__get_mp(),i,EA[i].__get_mp());
						  }
						  N[i]=0;
					  }else{
						  if(beverbose)
						  {
						  gmp_fprintf(stderr,"A[%d]=%Fe\t EA[%d]=%Fe\n",i,A[i].__get_mp(),i,EA[i].__get_mp());
						  }
						  N[i]=1;
						  l++;
					  }
				  }else
				  {
					  N[i]=0;
				  }
	  			}
			}

	  	  	//fprintf(stderr,"new range =%d\n",l);


			for(int i=0;i<k;i++)
			{
				Y[i]=0;
				for(int j=0;j<k;j++)
				{
					C[i][j]=0;
				}
			}

			sigm*=nsigm;
			nstars=crerat_matrix_sig(p,k,Nstars,avg_ef,massiv,C,Y,A,bitsize,sigm,&mag,N);

			*nstarss=nstars;

			mpf_class * AA = new mpf_class [l];
			mpf_class * YY = new mpf_class [l];
			mpf_class ** invCC = new mpf_class *[l];

			for(int i=0;i<l;i++)
			{
				invCC[i] = new mpf_class [l];
			}

			//cut range matrix and vector
			int li=0;
			for(int i=0;i<k;++i)
			{
				if(N[i]>0)
				{
				int lj=0;
				for(int j=0;j<k;++j)
				   {
					if(N[j]>0)
					{
					invCC[li][lj]=C[i][j];
					lj++;
					}
					//gmp_fprintf(stderr,"invC[%d][%d]=%Fe\n",i,j,invC[li][lj].__get_mp());
				   }
				YY[li]=Y[i];
				AA[li]=0;
				li++;
				}
			}

			inversion(invCC,l);

			multiplication(AA,invCC,YY,l);

			for(int i=0,j=0;i<k;++i)
			{
				if(N[i]>0)
				{
					invC[i][i]=invCC[j][j];
					A[i]=AA[j];
					j++;
				}else
				{
					A[i]=0.0;
					invC[i][i]=0.0;
				}
				//gmp_fprintf(stderr,"invC[%d][%d]=%Fe\n",i,i,invC[i][i].__get_mp());
			}


			calc_sigma(&sigm,&weightN,A,massiv,k,Nstars,avg_ef,p,bitsize,N);

			double sig=sqrt(mpf_get_d(sigm.__get_mp())/(weightN-l));

			  	if(beverbose)
			  	{
				fprintf(stderr,"Nstars=%ld weight=%lf sigma= %lf k=%d avg_mag=%6.3f\n",Nstars,weightN,sig,ici,mag);
			  	}


		/*	if(beverbose)
			{
				for(int i=0;i<k;++i)
				{

				  EA[i]=sqrt(sigm*invC[i][i]);
				  if( abs(A[i]/EA[i]) <= student_test)
				  {
					  gmp_fprintf(stderr,"A[%d]=%Fe\t EA[%d]=%Fe\n",i,A[i].__get_mp(),i,EA[i].__get_mp());
					  N[i]=0;
				  }else{
					  N[i]=1;
				  }
				}
			}
*/
			  	delete [] YY;
			  	delete [] AA;
			  	for(int i=0;i<l;++i)
			  	{
			  	delete [] invCC[i];
			  	}

			  	kk=l;
			  	ici++;
		  	}

	  	}


	  	*sigma=sigm;

	  	for(int i=0;i<k;++i)
	  	{
			delete [] C[i];
			delete [] invC[i];
	  	}
	  	delete [] Y;

	  	*magbin=mag;
return ici;

}

/*

*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static void show_usage( FILE * output, int argc, char * argv[] )
{
  fprintf(output, "Liniar regresion UTILITY\n");
  fprintf(output, "USAGE:\n");
  fprintf(output, "  %s OPTIONS FILE1\n", basename(argv[0]));
  fprintf(output, "  %s OPTIONS FILE2\n", basename(argv[1]));
  fprintf(output, "OPTIONS:\n");
  fprintf(output, "	x=		<integer>	one-based colum number of x column in file\n");
  fprintf(output, "	y=		<integer>	one-based colum number of y column in file, if don`t used, y=0\n");
  fprintf(output, "	mag=	<integer>	one-based colum number of mag column in file, if don`t used, mag=0\n");
  fprintf(output, "	f=		<integer>	one-based colum number of F column in file\n");
  fprintf(output, "	ef=		<integer>	one-based colum number of error F column in file, if used, else ef=0\n");
  fprintf(output, "	f1=		<integer>	one-based colum number of F1 column in file, if used\n");
  fprintf(output, "	ef1=		<integer>	one-based colum number of error F1 column in file, if used\n");
  fprintf(output, "	cx=		<integer>	colum number of corection x column in file2 if don`t used, if used\n");
  fprintf(output, "	cy=		<integer>	colum number of corection y column in file2 if  don`t used, y=0\n");
  fprintf(output, "	mag2=	<integer>	one-based colum number of mag column in file2, if don`t used, mag=0\n");
  fprintf(output, "	p=		<integer>	used polynomial degree, default p=1\n");
  fprintf(output, "	bitsize=<integer>	max bit size of float, default 256\n");
  fprintf(output, "	model=	<integer>	if used special reduction model, for ksi 1, for eta 2, don`t used default 0\n");
  fprintf(output, "	sigm=		<integer>	count simga for exclude, default 3sig\n");
  fprintf(output, "	ci=		<integer>	max count iteration, default 5\n");
  fprintf(output, "	st=		<double>	Students test, default 1\n");
  fprintf(output, "	-r					Print residuale in residual file\n");
  fprintf(output, "	-d					Print df in file2\n");
  fprintf(output, "	-c					Print coefficients\n");
  fprintf(output, "	-m					Print avg magnitude\n");
  fprintf(output, "	-n					Print count objects\n");
  fprintf(output, "	-e					Print coefficients error\n");
  fprintf(output, "	-v					Print some diagnostic messages to stderr (verbose mode)\n");
  fprintf(output, "	\n");

  UNUSED(argc);
}


int main(int argc, char *argv[])
{
	/* TODO: Fix this bug */

	const char * fname_in = {NULL};
	const char * fname_out = {NULL};

	FILE * fp; //FILE

	compression_t compression [2]={ compression_unknown, compression_unknown};


	int x=-1;		//	one-based colum number of x column in file
	int y=0;		//	one-based colum number of y column in file, if don`t used, y=0\n
	int mag=0;		//	one-based colum number of mag column in file, if don`t used, mag=0\n

	int cx=-1;		//	one-based colum number of x column in file 2
	int cy=-1;		//	one-based colum number of y column in file 2, if don`t used, y=0\n
	int mag2=0;		//	one-based colum number of mag column in file2, if don`t used, mag2=0\n


	int f=-1;		//	one-based colum number of F column in file
	int ef=0;		//	one-based colum number of error F column in file, if used, else ef=0\n
	int f1=0;		//	one-based colum number of F1 column in file
	int ef1=0;		//	one-based colum number of error F1 column in file, if used, else ef=0\n

	int p = 1; 		//	used polynomial degree, default p=1
	int nsigm = 3; 	//	count simga for exclude, default 3sig
	int ci = 5; 		//	max count iteration, default 5
	int bitsize = 256; 		//	max count iteration, default 256
	double st = 1; 		//	Students test, default 1
	int model = 0; 		//	Used, special model

	int print_coef 	= 0; 		//	Print coefficients
	int print_residual 	= 0; 		//	Print residual
	int print_error	= 0; 		//	Print coefficients error
	int beverbose 	= 0; 		//	Print some diagnostic messages to stderr (verbose mode)
	int usedfiles	= 0;		// Used File
	int print_df	= 0;		// Used File
	int avgmag		= 0;		//	Print avg mag
	int countobj	= 0;		//	Print count



	double avg_mag_bin=0;
	//double max_x = 0,min_x=1000000,max_y=0,min_y=1000000;

	/* parse command line */

	int i;
	for ( i = 1; i < argc; ++i )
	{
		 if ( strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-help") == 0 ) {
			  show_usage(stdout, argc, argv);
			  return 0;
			}

		/* read used columns from file1*/
		if ( strncmp(argv[i], "x=", 2) == 0 )
			{
			  if ( sscanf(argv[i] + 2, "%d", &x) != 1 || x < 1 )
			  {
				fprintf(stderr, "Invalid value of %s\n", argv[i]);
				return 1;
			  }
			}
		else if ( strncmp(argv[i], "y=", 2) == 0 )
		{
		  if ( sscanf(argv[i] + 2, "%d", &y) != 1)
		  {
			fprintf(stderr, "Invalid value of %s\n", argv[i]);
			return 1;
		  }
		}
		else if ( strncmp(argv[i], "mag=", 4) == 0 )
		{
		  if ( sscanf(argv[i] + 4, "%d", &mag) != 1)
		  {
			fprintf(stderr, "Invalid value of %s\n", argv[i]);
			return 1;
		  }
		}
		else if ( strncmp(argv[i], "f=", 2) == 0 )
		{
		  if ( sscanf(argv[i] + 2, "%d", &f) != 1 || f< 1 )
		  {
			fprintf(stderr, "Invalid value of %s\n", argv[i]);
			return 1;
		  }
		}
		else if ( strncmp(argv[i], "ef=", 3) == 0 )
		{
		  if ( sscanf(argv[i] + 3, "%d", &ef) != 1 )
		  {
			fprintf(stderr, "Invalid value of %s\n", argv[i]);
			return 1;
		  }
		}
		else if ( strncmp(argv[i], "f1=", 3) == 0 )
		{
		  if ( sscanf(argv[i] + 3, "%d", &f1) != 1 )
		  {
			fprintf(stderr, "Invalid value of %s\n", argv[i]);
			return 1;
		  }
		}
		else if ( strncmp(argv[i], "ef1=", 4) == 0 )
		{
		  if ( sscanf(argv[i] + 4, "%d", &ef1) != 1 )
		  {
			fprintf(stderr, "Invalid value of %s\n", argv[i]);
			return 1;
		  }
		}
		else if ( strncmp(argv[i], "cx=", 3) == 0 )
		{
		  if ( sscanf(argv[i] + 3, "%d", &cx) != 1 )
		  {
			fprintf(stderr, "Invalid value of %s\n", argv[i]);
			return 1;
		  }
		}

		else if ( strncmp(argv[i], "cy=", 3) == 0 )
		{
		  if ( sscanf(argv[i] + 3, "%d", &cy) != 1 )
		  {
			fprintf(stderr, "Invalid value of %s\n", argv[i]);
			return 1;
		  }
		}

		else if ( strncmp(argv[i], "mag2=", 5) == 0 )
		{
		  if ( sscanf(argv[i] + 5, "%d", &mag2) != 1)
		  {
			fprintf(stderr, "Invalid value of %s\n", argv[i]);
			return 1;
		  }
		}

		 /* read parameters for fitting */
		else if ( strncmp(argv[i], "p=", 2) == 0 )
		{
		  if ( sscanf(argv[i] + 2, "%d", &p) != 1 || p < 1 )
		  {
			fprintf(stderr, "Invalid value of %s\n", argv[i]);
			return 1;
		  }
		}
		else if ( strncmp(argv[i], "sigm=", 5) == 0 )
			{
			  if ( sscanf(argv[i] + 5, "%d", &nsigm) != 1 || nsigm < 1 )
			  {
				fprintf(stderr, "Invalid value of %s\n", argv[i]);
				return 1;
			  }
			}
		else if ( strncmp(argv[i], "ci=", 3) == 0 )
			{
			  if ( sscanf(argv[i] + 3, "%d", &ci) != 1 || ci < 1 )
			  {
				fprintf(stderr, "Invalid value of %s\n", argv[i]);
				return 1;
			  }
			}
		else if ( strncmp(argv[i], "st=", 3) == 0 )
			{
			  if ( sscanf(argv[i] + 3, "%lf", &st) != 1 || st < 1 )
			  {
				fprintf(stderr, "Invalid value of %s\n", argv[i]);
				return 1;
			  }
			}
		else if ( strncmp(argv[i], "bitsize=", 8) == 0 )
			{
			  if ( sscanf(argv[i] + 8, "%d", &bitsize) != 1 || bitsize < 1 )
			  {
				fprintf(stderr, "Invalid value of %s\n", argv[i]);
				return 1;
			  }
			}
		else if ( strncmp(argv[i], "model=", 6) == 0 )
			{
			  if ( sscanf(argv[i] + 6, "%d", &model) != 1 || model < 1 )
			  {
				fprintf(stderr, "Invalid value of %s\n", argv[i]);
				return 1;
			  }
			}


		else if ( strcmp(argv[i], "-v") == 0 ) {
			  beverbose = 1;
			}
		else if ( strcmp(argv[i], "-r") == 0 ) {
			  print_residual = 1;
			}
		else if ( strcmp(argv[i], "-c") == 0 ) {
			  print_coef = 1;
			}
		else if ( strcmp(argv[i], "-e") == 0 ) {
			  print_error = 1;
			}
		else if ( strcmp(argv[i], "-m") == 0 ) {
			  avgmag = 1;
			}
		else if ( strcmp(argv[i], "-n") == 0 ) {
			  countobj = 1;
			}
		else if ( strcmp(argv[i], "-d") == 0 ) {
			  print_df = 1;
			}
		else if ( !fname_in ) {
			  fname_in = argv[i];
			  usedfiles = 1;
			}
		else if ( !fname_out ) {
			  fname_out = argv[i];
			  usedfiles = 2;
			}
		else
			{
			  fprintf(stderr, "Invalid argument %s. Try %s --help\n", argv[i], argv[0]);
			  //printf("load data from file %s",fname);
			  return 1;
			}


	  }

	/* check command line inputs */
	  if ( !fname_in) {
		fprintf(stderr,"Input file names expected\n");
		show_usage(stderr, argc, argv);
		return -1;
	  }
	  if ( !fname_out && usedfiles==2) {
	  		fprintf(stderr,"Input file names expected\n");
	  		show_usage(stderr, argc, argv);
	  		return -1;
	  	  }

	  if ( x < 1 ) {
		  fprintf(stderr,"x argument is mandatory\n");
		  show_usage(stderr, argc, argv);
		  return -1;
		}
	  if ( y < 0 ) {
			  fprintf(stderr,"y argument is mandatory\n");
			  show_usage(stderr, argc, argv);
			  return -1;
			}
	  if ( mag < 0 ) {
			  fprintf(stderr,"mag argument is mandatory\n");
			  show_usage(stderr, argc, argv);
			  return -1;
			}

	  if ( cx < 1 && usedfiles==2) {
		  fprintf(stderr,"cx argument is mandatory\n");
		  show_usage(stderr, argc, argv);
		  return -1;
		}
	  if ( cy < 1 && usedfiles==2) {
			  fprintf(stderr,"cy argument is mandatory\n");
			  show_usage(stderr, argc, argv);
			  return -1;
			}
	  if ( mag2 < 0 ) {
			  fprintf(stderr,"mag2 argument is mandatory\n");
			  show_usage(stderr, argc, argv);
			  return -1;
			}

	  if ( f < 1 ) {
		  fprintf(stderr,"f argument is mandatory\n");
		  show_usage(stderr, argc, argv);
		  return -1;
		}
	  if ( ef < 0 ) {
		  fprintf(stderr,"ef argument is mandatory\n");
		  show_usage(stderr, argc, argv);
		  return -1;
		}
	  if ( f1 < 0 ) {
		  fprintf(stderr,"f1 argument is mandatory\n");
		  show_usage(stderr, argc, argv);
		  return -1;
		}
	  if ( ef1 < 0 ) {
		  fprintf(stderr,"ef1 argument is mandatory\n");
		  show_usage(stderr, argc, argv);
		  return -1;
		}

	  if ( p < 1 || p > 100) {
		fprintf(stderr,"p argument is mandatory\n");
		show_usage(stderr, argc, argv);
		return -1;
	  }

	  if ( nsigm < 1 || nsigm > 100) {
		fprintf(stderr,"sigm argument is mandatory\n");
		show_usage(stderr, argc, argv);
		return -1;
	  }

	  if ( ci < 1 || ci> 100) {
		fprintf(stderr,"ci argument is mandatory\n");
		show_usage(stderr, argc, argv);
		return -1;
	  }

	  if ( beverbose ) {
		  if(usedfiles==2)
		  {
		  fprintf(stderr,"Used file1 '%s' file2 '%s' and argumernts: x=%d\ty=%d\tf=%d\tef=%d\tp=%d\tsigm=%d\tci=%d\n",fname_in,fname_out,x,y,f,ef,p,nsigm,ci);
		  }else{
			  fprintf(stderr,"Used file '%s' and argumernts: x=%d\ty=%d\tf=%d\tef=%d\tf1=%d\tef1=%d\tp=%d\tsigm=%d\tci=%d\n",fname_in,x,y,f,ef,f1,ef1,p,nsigm,ci);
		  }
	  }

	  //printf("load data from file %s",fname);




	if (usedfiles>0)
	{
		/* check if input files are readable */
		    if ( access(fname_in, R_OK) != 0 ) {
		      fprintf(stderr, "Can't read %s: %s\n", fname_in, strerror(errno));
		      return -1;
		    }


		  if ( beverbose ) {
			fprintf(stderr,"loading %s....\n", fname_in);
		  }

		  if ( !(fp = open_file(fname_in, &compression[0])) ) {
			fprintf(stderr, "Can't read '%s': %s\n", fname_in, strerror(errno));
			return -1;
		  }

		  int k=0;	//count parameter fitting
		  if(y!=0)
		  {
			  for(int i=0;i<=p;++i)
			  {
				  k+=i+1;
			  }
		  }else
		  {k=p+1;}

		  //printf("k=%d\n",k);

		  obj_t *massiv = new obj_t[1000000];
		  obj_t *massiv1 = new obj_t[1000000];


		  /*
		  nstars =load_objects_xy(fp,x,y,f,ef,k,p,&C,&Y);
		  if (nstars==0 ) {
			fprintf(stderr, "Can't load %s\n", fname);
			return -1;
		  }
		  */

		  double avg_ef=0;
		  double avg_ef1=0;
		  long int Nstars=0;
		  if(f1>0)
		  {
			  long int temp_nstars=load_objects_xym(fp,x,y,mag,f,ef,&avg_ef,massiv);
			  fseek(fp,0,SEEK_SET);
			  Nstars=load_objects_xym(fp,x,y,mag,f1,ef1,&avg_ef1,massiv1);
			  if(temp_nstars!=Nstars)
			  {
				  fprintf(stderr, "ERROR in load data from %s Fdata = %ld and F1data = %ld\n", fname_in,temp_nstars,Nstars);
			  }
			  if(beverbose)
			  {
				  fprintf(stderr, "avg_ef=%lf  avg_ef1=%lf load %ld stars\n",avg_ef,avg_ef1,Nstars);
			  }

		  }else{
		  Nstars=load_objects_xym(fp,x,y,mag,f,ef,&avg_ef,massiv);
		  if(beverbose)
			  {
				  fprintf(stderr, "avg_ef=%lf\n",avg_ef);
			  }
		  }

		  close_file( fp, compression[0] );

		  if ( Nstars == 0 ) {
			fprintf(stderr, "Can't load data from %s \n", fname_in);
			return -1;
		  }

	  mpf_class *A 		= new mpf_class [k];
	  mpf_class *EA 	= new mpf_class [k];
	  mpf_class *A1 	= new mpf_class [k];
  	  mpf_class *EA1 	= new mpf_class [k];


	   for (int i = 0; i < k; i++)
	   {
	   A[i]=0;
	   EA[i]=0;
	   A1[i]=0;
	   EA1[i]=0;
	   }


  	mpf_class sigmax=1000000000000000;
  	mpf_class sigmay=1000000000000000;

  	int count=0;

  	if(f1>0)
	   {
		  	calc_polifit(&sigmax,nsigm,avg_ef,k,p,ci,bitsize,Nstars,massiv,A,EA,beverbose,&avg_mag_bin,&count,st,model);
		  	calc_polifit(&sigmay,nsigm,avg_ef1,k,p,ci,bitsize,Nstars,massiv1,A1,EA1,beverbose,&avg_mag_bin,&count,st,model);

		  	if(print_residual)
			{
		  		char fnamef[256];
		  		char fnamef1[256];
		  		sprintf(fnamef,"residual_f_%s.txt",fname_in);
		  		sprintf(fnamef1,"residual_f1_%s.txt",fname_in);
		  		FILE * frf=fopen(fnamef,"w");
				FILE * frf1=fopen(fnamef1,"w");
				pprint_residual(sigmax,nsigm,avg_ef,k,p,bitsize,Nstars,massiv,A,frf);
				pprint_residual(sigmay,nsigm,avg_ef1,k,p,bitsize,Nstars,massiv1,A1,frf1);
				fclose(frf);
				fclose(frf1);
			}
	   }else
	   {
		  	calc_polifit(&sigmax,nsigm,avg_ef,k,p,ci,bitsize,Nstars,massiv,A,EA,beverbose,&avg_mag_bin,&count,st,model);

		  	if(print_residual)
			{
		  	char *fnamef={NULL};
		  	sprintf(fnamef,"residual_f_%s.txt",fname_in);
			FILE * frf=fopen(fnamef,"w");
			pprint_residual(sigmax,nsigm,avg_ef,k,p,bitsize,Nstars,massiv,A,frf);
			fclose(frf);
			}
	   }

  	delete [] massiv;
  	delete [] massiv1;

  	if(usedfiles==2)
  	{

  		FILE *fr;
  		if ( !(fr = open_file(fname_out, &compression[1])) ) {
			fprintf(stderr, "Can't read '%s': %s\n", fname_out, strerror(errno));
			return -1;
		  }

  		int ncor=0;
  		if(f1>0)
  		{
		ncor=correction_objects(fr,cx,cy,f,f1,k,p,A,A1,bitsize,print_df);

  		}else{
  			ncor=correction_objects(fr,cx,cy,f,k,p,A,bitsize,print_df);
  		}

  		fprintf(stderr,"Correction %d objects in file %s\n", ncor,fname_out);
		close_file( fr, compression[1] );

  	}else{
  	if(f1>0)
  	{
  		if(avgmag==1)
  		{
  		printf("%10.3f\t",avg_mag_bin);
  		}
  		if(print_coef==1)
			{
				for(int i=0;i<k;++i)
				  {
				  //sigma=sqrt(sigma);
				  gmp_printf("%Fe\t",A[i].__get_mp());
				  }
				printf("\n");
			}
		if(print_error==1)
		{
			for(int i=0;i<k;++i)
			  {
			  //sigma=sqrt(sigma);
			  gmp_printf("%Fe\t",EA[i].__get_mp());
			  }
			printf("\n");
		}
		if(print_coef==1)
		{
			for(int i=0;i<k;++i)
			  {
			  //sigma=sqrt(sigma);
			  gmp_printf("%Fe\t",A1[i].__get_mp());
			  }
			printf("\n");
		}
		if(print_error==1)
		{
			for(int i=0;i<k;++i)
			  {
			  //sigma=sqrt(sigma);
			  gmp_printf("%Fe\t",EA1[i].__get_mp());
			  }
			printf("\n");
		}
  	}else
  	{
  		if(avgmag==1)
		{
		printf("%10.3f\t",avg_mag_bin);
		}

  		if(print_coef==1)
		{
			for(int i=0;i<k;++i)
			  {
			  //sigma=sqrt(sigma);
			  gmp_printf("%Fe\t",A[i].__get_mp());
			  }
			if(countobj==1)
					{
						printf("%d\n",count);
					}else{	printf("\n");}
		}
		if(print_error==1)
		{
			for(int i=0;i<k;++i)
			  {
			  //sigma=sqrt(sigma);
			  gmp_printf("%Fe\t",EA[i].__get_mp());
			  }if(countobj==1)
				{
					printf("%d\n",count);
				}else{	printf("\n");}
		}

  	}
  	}

	if ( beverbose )
	{
		if(f1>0)
		{
		for(int i=0;i<k;++i)
		  {
		  //sigma=sqrt(sigma);
		  gmp_fprintf(stderr,"A[%d]=%Fe +-%Fe\t A1[%d]=%Fe +-%Fe\n",i,A[i].__get_mp(),EA[i].__get_mp(),i,A1[i].__get_mp(),EA1[i].__get_mp());
		  }
		printf("\n");
		}else{
			for(int i=0;i<k;++i)
			  {
			  //sigma=sqrt(sigma);
			  gmp_fprintf(stderr,"A[%d]=%Fe +-%Fe\n",i,A[i].__get_mp(),EA[i].__get_mp());
			  }

		}
	}
	delete [] A;
	delete [] A1;
	delete [] EA;
	delete [] EA1;
	}

	return 0;
}
