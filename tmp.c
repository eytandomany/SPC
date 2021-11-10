#include "SW.h"

#define S 10000

void main() {
  int i,k=0,indx[S];
  double a[S];
  for( i=0;i<S;i++ ) a[i]=cos(i);
  DSortIndex(S,a,indx);
  for( i=1;i<S;i++ )
     if( a[indx[i]]<=a[indx[i-1]] )
        k++;
  printf("missed %d out of %d\n",k,S);
}
