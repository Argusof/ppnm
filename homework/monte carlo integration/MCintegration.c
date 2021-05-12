#include<math.h>
#include<stdlib.h>
#include<assert.h>

#define RND ((double)rand()/RAND_MAX)

void randompts(int dim, double* lowerBound, double* upperBound, double* pts){
  for(int i = 0; i < dim; i++){
    pts[i] = lowerBound[i] + RND*(upperBound[i] - lowerBound[i]);
  }
}

void plainMC(int dim, double* lowerBound, double* upperBound, double f(double* pts), int numOfPts, double* result, double* error ){
  double volume = 1;
  for(int i = 0; i < dim; i++){
    volume *= upperBound[i] - lowerBound[i];
  }

    double sum = 0;
    double sumSquare = 0;
    double funcVal;
    double pts[dim];

    for(int i = 0; i < numOfPts; i++){
      randompts(dim, lowerBound, upperBound, pts);
      funcVal = f(pts);
      sum += funcVal;
      sumSquare += funcVal*funcVal;
    }

    double avr = sum/numOfPts;
    double var = sumSquare/numOfPts - avr*avr;
    *result = avr*volume;
    *error = sqrt(var/numOfPts)*volume;
}


double corput(int id, int base){
  double corputNum = 0;
  double coprimeBase = (double)1/base;
  while(id > 0){
    corputNum += (id%base)*coprimeBase;
    id /= base;
    coprimeBase /= base;
  }
  return corputNum;
}

void halton(int id, int dim, double* pts){
  int base[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79};
  int maptsDim = sizeof(base)/sizeof(int);
  assert(dim <= maptsDim );
  for(int i = 0; i < dim; i++){
    pts[i] = corput(id + 1, base[i]);
  }
}

void halton2(int id, int dim, double* pts){
  int base[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83};
  int maptsDim = sizeof(base)/sizeof(int);
  assert(dim <= maptsDim );
  for(int i = 0; i < dim; i++){
    pts[i] = corput(id + 1, base[i]);
  }
}

void HaltonCorputRandompts(int id, int dim, double* lowerBound, double* upperBound, double* pts){
  halton(id, dim, pts);
  for(int i = 0; i < dim; i++){
    pts[i] = lowerBound[i] + pts[i]*(upperBound[i] - lowerBound[i]);
  }
}

void HaltonCorputRandompts2(int id, int dim, double* lowerBound, double* upperBound, double* pts){
  halton2(id, dim, pts);
  for(int i = 0; i < dim; i++){
    pts[i] = lowerBound[i] + pts[i]*(upperBound[i] - lowerBound[i]);
  }
}

void HaltonCorputMC(int dim, double* lowerBound, double* upperBound, double f(double* pts), int numOfPts, double* result, double* error ){
  double volume = 1;
  for(int i = 0; i < dim; i++){
    volume *= upperBound[i] - lowerBound[i];
  }

    double sum = 0;
    double sum2 = 0;
    double sumSquare = 0;
    double sumSquare2 = 0;
    double funcVal;
    double funcVal2;
    double pts[dim];
    double pts2[dim];

    for(int i = 0; i < numOfPts/2; i++){
      HaltonCorputRandompts(i, dim, lowerBound, upperBound, pts);
      HaltonCorputRandompts2(i, dim, lowerBound, upperBound, pts2);
      funcVal = f(pts);
      funcVal2 = f(pts2);
      if(!isinf(funcVal) && !isinf(funcVal2)){
        sum += funcVal;
        sum2 += funcVal2;
      }
    }

    *result = (sum+sum2)/numOfPts*volume;
    *error = fabs(sum-sum2)/numOfPts*volume;
}


double stratMC(int dim, double f(double* pts), double* lowerBound, double* upperBound, double absAcc, double relAcc, int numOfRecalls, double meanRecall){

  int numOfPts = 16*dim;
  double volume = 1;
  for(int k = 0; k < dim; k++){
     volume *= upperBound[k] - lowerBound[k];
  }
  int numOfPtsLeft[dim];
  int numOfPtsRight[dim];
  double pts[dim];
  double meanleft[dim];
  double meanright[dim];
  double mean = 0;
  for(int k = 0; k < dim; k++){
    meanleft[k]      = 0;
    meanright[k]     = 0;
    numOfPtsLeft[k]  = 0;
    numOfPtsRight[k] = 0;
  }
  for(int i = 0; i < numOfPts; i++){
      //pts[k] = lowerBound[k] + RND*(upperBound[k] - lowerBound[k]);
      randompts(dim, lowerBound, upperBound, pts);
      double funcVals = f(pts);

      for(int k = 0; k < dim; k++){
        if(pts[k] > (lowerBound[k] + upperBound[k])/2){
          numOfPtsRight[k]++;
          meanright[k] += funcVals;
        }
        else{
          numOfPtsLeft[k]++;
          meanleft[k] += funcVals;
        }
      }
      mean += funcVals;
    }
    mean /= numOfPts;

    for(int k = 0; k < dim; k++){
      meanleft[k]  /= numOfPtsLeft[k];
      meanright[k] /= numOfPtsRight[k];
    }
    int kdiv = 0;
    double maxvar = 0;
    for(int k = 0; k < dim; k++){
      double var = fabs(meanright[k] - meanleft[k]);
      if(var > maxvar){
        maxvar = var;
        kdiv = k;
      }
    }

    double result = (mean*numOfPts + meanRecall*numOfRecalls)/(numOfPts + numOfRecalls)*volume;
    double error = fabs(meanRecall - mean)*volume;
    double tolerance = absAcc + fabs(result)*relAcc;
    if(error < tolerance){
      return result;

    }
    double lowerBound2[dim];
    double upperBound2[dim];
    for(int k = 0; k < dim; k++){
      lowerBound2[k] = lowerBound[k];
      upperBound2[k] = upperBound[k];
    }
    lowerBound2[kdiv] = (lowerBound[kdiv] + upperBound[kdiv])/2;
    upperBound2[kdiv] = (lowerBound[kdiv] + upperBound[kdiv])/2;

    double resultLeft  = stratMC(dim, f, lowerBound, upperBound2, absAcc/sqrt(2), relAcc, numOfPtsLeft[kdiv], meanleft[kdiv]);
    double resultRight = stratMC(dim, f, lowerBound2, upperBound, absAcc/sqrt(2), relAcc, numOfPtsRight[kdiv], meanright[kdiv]);

    return resultLeft + resultRight;
  }



double stratamc(int dim, double f(int dim, double* x), double* a, double* b,double acc, double eps, int n_reuse, double mean_reuse){

int N=16*dim;
double V=1;
for (int i=0;i<dim;i++){
	V*=b[i]-a[i];
}

int n_l[dim];
int n_r[dim];

double x[dim];
double mean_l[dim];
double mean_r[dim];
double mean=0;

unsigned int rand_seed=3;


for (int i=0;i<N;i++){
	for(int j=0;j<dim;j++){
		x[j]=a[j]+rand_r(&rand_seed)*1.0/RAND_MAX*(b[j]-a[j]);
	}

	double fx=f(dim,x);
	for (int j=0;j<dim;j++){
		if(x[j]>(a[j]+b[j])/2){
			n_r[j]++;
			mean_r[j]+=fx;
		} else {
			n_l[j]++;
			mean_l[j]+=fx;
		}
	}

	mean+=fx;
}
mean/=N;

for(int i=0;i<dim;i++){
	mean_l[i]/=n_l[i];
	mean_r[i]/=n_r[i];
}

int idiv=0;
double maxvar=0;

for(int i =0;i<dim;i++){
	double var=fabs(mean_r[i]-mean_l[i]);
	if(var>maxvar){
		maxvar=var;
		idiv=i;
	}
}


double integ=(mean*N+mean_reuse*n_reuse)/(N+n_reuse)*V;
double err=fabs(mean_reuse-mean)*V;
double tol=acc+fabs(integ)*eps;
if (err<tol){
	return integ;
}

double a2[dim];
double b2[dim];

for (int i=0;i<dim;i++){
	a2[i]=a[i];
	b2[i]=b[i];
}

a2[idiv]=(a[idiv]+b[idiv])/2;
b2[idiv]=(a[idiv]+b[idiv])/2;
double integ_l=stratamc(dim,f,a,b2,acc/sqrt(2.0),eps,n_l[idiv],mean_l[idiv]);
double integ_r=stratamc(dim,f,a2,b,acc/sqrt(2.0),eps,n_r[idiv],mean_r[idiv]);
return integ_l+integ_r;

}
