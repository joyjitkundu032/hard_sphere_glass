/*CHANGE THE FUNCTION "read_input" if you change the dimensions D. Here Rmin, rd correspond to the minimum and maximum diameter of the spheres respectively.*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <time.h>
#include "ran2.c"

#define eps 0.0000000001
#define D 6
#define Ravg 1.0
#define f 0.027
#define tr (f*Ravg)
#define swapeq 0.20
#define TD 1000
#define swapprob 0.00
#define IFLAG 1
#define dsigmatol 0.090
#define pi (22.0/7.0)

int N,Teq,nsteps,GAP,INIT,MNNEI;
double boxsize,cellsize,boxl,boxl2,rd,Rmin,vf,drneimax,drneimax2,Rskin,Rskin_in;
char outfile4[200],readfile[200],outfile5[200],outfile6[200];
long int seed=485620;
double *RD,*dispsum[D],dCOM[D];
int **NNList,*Ncross[D];

FILE *fpw;

/*this program returns the nearest integer*/
double anint(double x)
{
	double d;
	if(x>=0.50)
		d=1.00;
	else
	{
		if(x<=(-0.50))
			d=-1.00;
		else
			d=0.00;
	}
	return d;
} /*end of the anint body*/

void read_input(double r[D][N], double VF)
{
	FILE *fpr;
	int i,k,time;
	time=INIT;
	int c1,c2,c3,c4,c5,c6;
	double c7,c8,c9,c10,c11,c12;
	i=0;
	if(IFLAG == 1)
	{
		sprintf(readfile,"./../contract_out/fconfig%dD_N%d_BS%1.4lfRmax%1.3lfRmin%1.3lf_vf%1.4lf.dat",D,N,boxsize,rd,Rmin,VF);
		fpr=fopen(readfile,"r");
		while(fscanf(fpr,"%lf%lf%lf%lf%lf%lf%lf",&r[0][i],&r[1][i],&r[2][i],&r[3][i],&r[4][i],&r[5][i],&RD[i])!=EOF)
			i++;
	}
	else if(IFLAG == 2)
	{
                sprintf(readfile,"rho_%1.4lfSWP_0.20eq_seed/config%dD_N%dBS%1.4lfRmax%1.3lfRmin%1.3lf_vf%1.4lfRskin_%1.2lft_%d_10.dat",vf,D,N,boxsize,rd,Rmin,vf,Rskin_in,time);
		fpr=fopen(readfile,"r");
		while(fscanf(fpr,"%lf%lf%lf%lf%lf%lf%d%d%d%d%d%d%lf%lf%lf%lf%lf%lf%lf",&r[0][i],&r[1][i],&r[2][i],&r[3][i],&r[4][i],&r[5][i],&c1,&c2,&c3,&c4,&c5,&c6,&c7,&c8,&c9,&c10,&c11,&c12,&RD[i])!=EOF)
			i++;
	}
	fclose(fpr);
}

void next(int v[D],int m)
{
        int ip;
        ip=0;
        while(v[ip] < m)
        {
                v[ip]=v[ip]+1;
                if(v[ip] < m)
                        break;
                else if(v[ip] > (m-1))
                {
                        v[ip]=0;
                        ip=ip+1;
                }
        }
}

void create_neighbour_list(double r[D][N])
{
	int k,j,i,nnei,nnlistbeg;
        double rdiff[D],rdiffsq,rdiffsqrt,sigma;

	for(i=0;i<N;i++)
	{
		nnei=0;
		sigma=RD[i];
		for(j=0;j<N;j++)
		{
			if(j != i)
			{
				rdiffsq=0.0;
				for(k=0;k<D;k++)
				{
					rdiff[k]=r[k][i]-r[k][j];
					rdiff[k]=rdiff[k]-(boxl)*anint(rdiff[k]/(boxl));
					rdiffsq=rdiffsq+rdiff[k]*rdiff[k];
				}
				rdiffsqrt=sqrt(rdiffsq);
				if(rdiffsqrt < cellsize)
				{
					nnei++;
					NNList[i][nnei]=j;

				}
			}
		}
		NNList[i][0]=nnei;
	}
}

int ipow(int base, int exp)
{
	int result = 1;
	while (exp)
	{
		if (exp & 1)
			result *= base;
		exp >>= 1;
		base *= base;
	}
	return result;
}

double correct_pbc(double rr)
{
        double rr_n;
        if(rr >= boxl2)
                rr_n=rr-boxl;
        else if(rr < -boxl2)
                rr_n=rr+boxl;
        else if(rr > (boxl2-eps))
                rr_n=-boxl2;
        else
                rr_n=rr;
        return rr_n;
}

int check_overlap(int j, double posi[D], double r[D][N], double sigma)
{
	int k,tmpj,nb;
	double rdiff[D],rdiffsq,rdiffsqrt;
	
	for(nb=1;nb<=NNList[j][0];nb++)
	{
		tmpj=NNList[j][nb];
		rdiffsq=0.0;
		for(k=0;k<D;k++)
		{
			rdiff[k]=posi[k]-r[k][tmpj];
			rdiff[k]=rdiff[k]-(boxl)*anint(rdiff[k]/(boxl));
			rdiffsq=rdiffsq+rdiff[k]*rdiff[k];
		}
		rdiffsqrt=sqrt(rdiffsq);
		if(rdiffsqrt <((sigma+RD[tmpj])/2.0+eps))
			return 1;
	}	
	return 0;
}

int check_overlap_swap(int j, int j2, double posi[D], double r[D][N], double sigma)
{
        int k,tmpj,nb;
        double rdiff[D],rdiffsq,rdiffsqrt;

        for(nb=1;nb<=NNList[j][0];nb++)
        {
                tmpj=NNList[j][nb];
		//if(ic > 511)
			//printf("NN[%d][%d]=%d tmpj=%d\n",ic,nb,NN[ic][nb],tmpj);
		if(tmpj != j && tmpj != j2)
		{
			rdiffsq=0.0;
			for(k=0;k<D;k++)
			{
				rdiff[k]=posi[k]-r[k][tmpj];
				rdiff[k]=rdiff[k]-(boxl)*anint(rdiff[k]/(boxl));
				rdiffsq=rdiffsq+rdiff[k]*rdiff[k];
			}
			rdiffsqrt=sqrt(rdiffsq);
			if(rdiffsqrt <((sigma+RD[tmpj])/2.0+eps))
				return 1;
		}      
	}
        return 0;
}

int check_overlap_N(int j, double posi[D], double r[D][N], double sigma)
{
        int k,nb;
        double rdiff[D],rdiffsq,rdiffsqrt;

        for(nb=0;nb<N;nb++)
        {
		if(nb != j)
		{
			rdiffsq=0.0;
			for(k=0;k<D;k++)
			{
				rdiff[k]=posi[k]-r[k][nb];
				rdiff[k]=rdiff[k]-(boxl)*anint(rdiff[k]/(boxl));
				rdiffsq=rdiffsq+rdiff[k]*rdiff[k];
			}
			rdiffsqrt=sqrt(rdiffsq);
			if(rdiffsqrt <((sigma+RD[nb])/2.0+eps))
				return 1;
		}
	}
        return 0;
}

void find_max_displacement(double dr)
{
	if(dr > drneimax)
	{
		drneimax2=drneimax;
		drneimax=dr;
	}
	else
		if(dr > drneimax2)
			drneimax2=dr;
}

void clear_displacement()
{
	int i,k;
	for(i=0;i<N;i++)
		for(k=0;k<D;k++)
			dispsum[k][i]=0.0;
	drneimax=0.0; drneimax2=0.0;
} 

void restore_dis(int i, double disp[D], double x1, double x2)
{
	int k;
	for(k=0;k<D;k++)
		dispsum[k][i]=dispsum[k][i]-disp[k];
	drneimax2=x2; drneimax=x1;
}

void translate(double r[D][N])
{
        int k,param,param1,j,i;
        double disp[D],tmpr[D],input_r,sigma,drneisq,drnei,tmp1,tmp2,cross[D],ur,true_r;
        j=(int)(ran2(&seed)*N);

        for(k=0;k<D;k++)
        {
                disp[k]=2.0*tr*ran2(&seed)-tr;
                input_r=r[k][j]+disp[k];
		cross[k]=input_r;
                tmpr[k]=correct_pbc(input_r);
        }
        sigma=RD[j];
        param=check_overlap(j,tmpr,r,sigma);
        if(param != 0)
                return;
        drneisq=0.0;
        for(k=0;k<D;k++)
        {
                dispsum[k][j]=dispsum[k][j]+disp[k];
                drneisq=drneisq+pow(dispsum[k][j],2.0);
        }
        drnei=sqrt(drneisq);
        tmp1=drneimax; tmp2=drneimax2;
        find_max_displacement(drnei);
        if((drneimax+drneimax2) > (Rskin-eps))
        {
                param1=check_overlap_N(j,tmpr,r,sigma);
                if(param1 == 0)
                {
                        for(k=0;k<D;k++)
                        {
                                r[k][j]=tmpr[k];
                                dCOM[k]=dCOM[k]+disp[k];
                                if((cross[k]) > boxl2)
                                        Ncross[k][j]=Ncross[k][j]+1;
                                else if ((cross[k]) < -boxl2)
                                        Ncross[k][j]=Ncross[k][j]-1;
                        }
                        create_neighbour_list(r);
                        clear_displacement();
		}
                else
                        restore_dis(j,disp,tmp1,tmp2);
        }
        else
        {
                for(k=0;k<D;k++)
                {
                        r[k][j]=tmpr[k];
                        dCOM[k]=dCOM[k]+disp[k];
                        if((cross[k]) > boxl2)
                                Ncross[k][j]=Ncross[k][j]+1;
                        else if ((cross[k]) < -boxl2)
                                Ncross[k][j]=Ncross[k][j]-1;
                }
        }
}

void swap(double r[D][N])
{
	int i1,i2,param1,param2,k;
	double sigma1,sigma2,dsigma,tmpr1[D],tmpr2[D];

	dsigma=1.0;		
	while(dsigma > dsigmatol)
	{
		i1=(int)(ran2(&seed)*N);
		i2=(int)(ran2(&seed)*N);
		if(i1 != i2)
			dsigma=fabs(RD[i1]-RD[i2]);
	}

	for(k=0;k<D;k++)
	{
		tmpr1[k]=r[k][i1];
		tmpr2[k]=r[k][i2];
	}
	sigma1=RD[i1]; sigma2=RD[i2];
	param1=check_overlap_swap(i1,i2,tmpr1,r,sigma2);
	param2=check_overlap_swap(i2,i1,tmpr2,r,sigma1);
	if((param1+param2) == 0)
	{
		RD[i1]=sigma2; RD[i2]=sigma1;
	}
}

void evolve(double r[D][N], double pswp)
{
	int i;
	//printf("I AM HERE\n");
	for(i=0;i<N;i++)
	{
		if(ran2(&seed) < pswp)
			swap(r);
		else
			translate(r);
	}
}
double vol_measure(double diam)
{
        double vol,g,dd,ag;
        dd=1.0*D;
        ag=D/2.0+1.0;
        vol=pow(pi,dd/2.0)*pow(diam,dd)/tgamma(ag);
        return vol;
}

double measure_packing_frac()
{
        int i;
        double vol,est;
        vol=0.0;
        for(i=0;i<N;i++)
                vol=vol+1.0*vol_measure(RD[i]);
        vol=vol/pow(boxl,1.0*D)/ipow(2,D);
        printf("packing fraction = %lf\n",vol);
        return vol;
}

void print_configuration(double r[D][N], int time)
{
        int i,k;
        double input_r,tmpr[D],rr[D];
        sprintf(outfile4,"./rho_%1.4lfSWP_%1.2lfeq/config%dD_N%dBS%1.4lfRmax%1.3lfRmin%1.3lf_vf%1.4lfRskin_%1.2lft_%d.dat",vf,swapprob,D,N,boxsize,rd,Rmin,vf,Rskin,time);
        fpw=fopen(outfile4,"w");
        for(i=0;i<N;i++)
        {
                for(k=0;k<D;k++)
                {
                        input_r=r[k][i]-1.0*dCOM[k]/N;
                        rr[k]=input_r;
                        tmpr[k]=correct_pbc(input_r);
                        fprintf(fpw,"%.16e\t",tmpr[k]);
                }
                for(k=0;k<D;k++)
                        fprintf(fpw,"%d\t",Ncross[k][i]);
                for(k=0;k<D;k++)
                        fprintf(fpw,"%.16e\t",rr[k]+boxl*Ncross[k][i]);
                fprintf(fpw,"%.16e\n",RD[i]);
        }
        fclose(fpw);
}

int main(void)
{
	FILE *fp,*fpt;
	int i,it;
	time_t  start_time = clock();

	scanf("%d",&N);
	scanf("%lf",&boxsize);
	scanf("%lf",&rd);
	scanf("%lf",&Rmin);
	scanf("%d",&nsteps);
	scanf("%d",&Teq);
	scanf("%lf",&vf);
	scanf("%d",&INIT);
	scanf("%lf",&Rskin);
	scanf("%lf",&Rskin_in);
	/*N=3400; boxsize=1.563903; rd=1.50241; Rmin=0.7494; 
	nsteps=1000; Teq=0; INIT=799200; vf=0.3620;*/
	boxl=boxsize; boxl2=((boxl)/2.0); GAP=(nsteps/TD);
	cellsize=(rd+Rskin);
	MNNEI=N;
	float time1 = (float) (clock() - start_time) / CLOCKS_PER_SEC;

	RD=(double *) malloc (N * sizeof(double));
	NNList=(int **) malloc (N * sizeof(int *));
	for(it=0;it<D;it++)
	{
		dispsum[it]=(double *) malloc (N * sizeof(double));
		Ncross[it]=(int *) malloc (N * sizeof(int));
	}
	for(it=0;it<N;it++)
		NNList[it]=(int *) malloc (MNNEI * sizeof(int));
	double r[D][N];
	sprintf(outfile5,"time_tracker_%dD_N%dBS%1.4lfRmax%1.3lfRmin%1.3lf_vf%1.4lfSWP%1.2lfRskin%1.2lf.dat",D,N,boxsize,rd,Rmin,vf,swapprob,Rskin);
	sprintf(outfile6,"clock_time_%dD_N%dBS%1.4lfRmax%1.3lfRmin%1.3lf_vf%1.4lfSWP%1.2lfRskin%1.2lf.dat",D,N,boxsize,rd,Rmin,vf,swapprob,Rskin);
	read_input(r,vf);
	
	measure_packing_frac();
	create_neighbour_list(r);

	for(it=0;it<Teq;it++)
	{
		evolve(r,swapeq);
		if(it % GAP == 0)
		{
			fp=fopen(outfile5,"w");
			fprintf(fp,"# it = %d\n",it);
			fclose(fp);	
		}	
	}	
	
	for(it=0;it<nsteps;it++)
	{
		evolve(r,swapprob);
		if(it % GAP == 0)
			print_configuration(r,it);
	}
	float time2 = (float) (clock() - start_time) / CLOCKS_PER_SEC;
        fpt=fopen(outfile6,"w");
        fprintf(fpt,"Time taken=%f seconds\n", time2-time1);
}
