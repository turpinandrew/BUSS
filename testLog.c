#include <stdio.h>
#include <stdlib.h>
#include <math.h>



int LogTable256[256];

inline double myLog2(double p) {
    int v = (int) ( (double)p * (double)(1<< 30));
    int tt, r;

printf("%d\n",v);
    if (tt = v >> 24) 
        r = 24 + LogTable256[tt];
    else if (tt = v >> 16) 
        r = 16 + LogTable256[tt];
    else if (tt = v >> 8) 
        r = 8 + LogTable256[tt];
    else 
        r = LogTable256[v];

    return (double)r - 30.0;
}//myLog2()


inline float fast_log2 (float val)
{
   int * const    exp_ptr = (int *)(&val);
   int            x = *exp_ptr;
   const int      log_2 = ((x >> 23) & 255) - 128;
   x &= ~(255 << 23);
   x += 127 << 23;
   *exp_ptr = x;

   val = ((-1.0f/3) * val + 2) * val - 2.0f/3;   // (1)

   return (val + log_2);
} 

int
main() {

    LogTable256[0] = LogTable256[1] = 0;
    for (int i = 2; i < 256; i++) 
        LogTable256[i] = 1 + LogTable256[i / 2];
    LogTable256[0] = -1; 

    for(int i = 1 ; i < 10 ; i++) {
        double p = (double)i / 10.0;
        printf("%10.8f %10.8f %10.8f \n", log2(p), fast_log2(p), myLog2(p));
    }
}

