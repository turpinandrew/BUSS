/*
** A C version of BUSS.
** Based on surface5_3.r, but without RLE and Memoization.
**
** The i'th value of the pdf/lf domain, x, is the base |D| representation
** of i, and the digits of that value are indexes into D.
** So if D = {1,17,31} and |L| = 4
**  i =  0  x = 0000  sample = { 1, 1, 1, 1}
**  i =  8  x = 0022  sample = { 1, 1,31,31}
**  i = 25  x = 0221  sample = { 1,31,31,17}  ( 2*3^2 + 2*3^1 + 2*3^0 )
**  i = 57  x = 2011  sample = {31, 1,17,17}  ( 2*3^3 + 0*3^2 + 1*3^1 + 1*3^0 )
**
** Andrew Turpin
** Wed  2 Mar 2011 13:14:19 EST
** Modified Wed  9 Mar 2011: changed to linked tree structure for memo.
**                           code has a horrible mix of chunks and trees.
**                           Really need to abstract TreeNode * out of loop in buss(). Done inefficiently.
** Modified Thu 10 Mar 2011: Remove pdfs from tree and regen on demand.
** Modified April 2012:      Added bimodal pdfs as priors
**                           Removed insistance on not presenting at same location twice in a row.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include "smalloc.h"
#include "ent.h"
#include "chunk.h"
#include "tree.h"

// #include "log.h"
// #define LOG(_p) log2(_p)
//#define LOG(_p) fast_log2(_p)

typedef float Prob; // in case I want double later...

typedef struct moments { 
    Prob mean, sd, H;
} MomentType; 

typedef struct fos { 
    float threshold, fp, fn, slope; 
} Fos;

#define min(_a, _b) ((_a) < (_b) ? (_a) : (_b))

    // prob of seeing _x with FoS _fos
#define PSI(_x, _fos) ((_fos).fp + (1-(_fos).fn-(_fos).fp)*(gsl_cdf_gaussian_Q((_x)-(_fos).threshold, (_fos).slope)))

   // name of prior file to read if priorFlag == 0
#define PRIOR_FILENAME "buss_prior.txt"

   // random numbers
gsl_rng *rng;

    // command line parameters
int numberOfTrials;
int numberOfLocations;
int numberOfDbValues;
int *dbValues;
Fos *fos;
int priorType;
int lfType;
char stopType;
float stopValue;
int maxDepth;
char rewriteTreeFlag;   // if true, rewrite tree at the end if it has changed
char verbose;   // flag

    // and derived form command line params
char (*stopFunction)(Prob *pdf, int depth);
int domain;   // |D|^|L| = numberOfDbValues^numberOfLocations
Prob *prior;  // Pr(domain is true surface sample)
Prob **lfLut; // lfLut[s, t] = Pr(yes to dbValues[s] | tt=dbValues[t])
char **digits; // digits[i][l] is the l'th digit (LH=0,RH=numberOfLocations-1) of the base |D| expansion of i \in[0..domain-1]
               // Assumes |D| does not exceed 127

    // For decision tree memoization 
Tree *tree;
char treeChanged;   // global flag to see if worth rewriting tree at the end.

#define TO_N_PTR(_i) ((TreeNode *)chunkGet(tree->nodes, _i))    // bit inefficient, but scary realloc!

Prob *workArea1; // working area for pdf*lf multiplication in chooseStimulus()
Prob *workArea2; // for creating nodes in makeNewNode()
Prob *workArea3; // for creating pdf when calling makePdf()

    // prototypes
char stopStdev(Prob *pdf, int depth)   ;
char stopEntropy(Prob *pdf, int depth) ;
char stopNumPres(Prob *pdf, int depth) ;

/*
** Return a LUT lf[stimIndex, trueThreshIndex] 
**    = Pr(yes to dbValues[stimIndex] | tt=dbValues[trueThreshIndex])
*/
void
makeLf(int lfType) {
    lfLut = (Prob **)smalloc(sizeof(Prob *) * numberOfDbValues);

    for(int s = 0 ; s < numberOfDbValues ; s++) {
        lfLut[s] = (Prob *)smalloc(sizeof(Prob) * numberOfDbValues);
        for(int tt = 0 ; tt < numberOfDbValues ; tt++) {
            if (lfType == 1) 
                lfLut[s][tt] = 0.03 + (1-0.03-0.03)*(gsl_cdf_gaussian_Q(dbValues[s]-dbValues[tt], min(1.0, exp(-0.098 * tt +
3.62))));
//printf("%2d %2d %f\n",dbValues[s],dbValues[tt], lfLut[s][tt]);
        }
    }
}//makeLf()

// Floor of PRIOR_4_FLOOR + line (t-L,f)->(t,H), line (t,H)->(t+M,f)
// ie box with skew triangle on top
#define PRIOR_4_L 10    // number of db steps to the lower end to linear over
#define PRIOR_4_M 3     // number of db steps to the upper end to linear over
#define PRIOR_4_FLOOR 0.00000001  // smallest prob for low end (prior to normalisation)
Prob 
prior4(char t1, char t2) {
    char M = (Prob)min(PRIOR_4_M, 40-t1);
    char L = (Prob)min(PRIOR_4_L, t1);
    Prob H = 2.0*(1 - 40.0*PRIOR_4_FLOOR)/(Prob)(L+M) + PRIOR_4_FLOOR;
//printf("%2d %2d %2d %2d %f %f\n",t1,t2,L,M,H);
    if ((t2 < t1 - L) || (t2 > t1 + M)) 
        return (Prob)PRIOR_4_FLOOR;
    else {
        if (t2 <= t1) {
            return H - (L>0 ? (H - PRIOR_4_FLOOR) * (Prob) (t1 - t2) / (Prob)L : (Prob)0.0);
        } else {
            return H - (M>0 ? (H - PRIOR_4_FLOOR) * (Prob) (t2 - t1) / (Prob)M : (Prob)0.0);
        }
    }
}//prior4()

// Floor of PRIOR_5_FLOOR + line (t-L,f)->(t,H), line (t,H)->(t+M,f)
// ie box with skew triangle on top
// triangle has height topLevel, and floor of pdf is 1 (prior to normalisation)
#define PRIOR_5_L       10 // number of db steps to the lower end to linear over
#define PRIOR_5_M        3 // number of db steps to the upper end to linear over
#define PRIOR_5_FLOOR    1 // smallest "prob" (prior to normalisation)
#define PRIOR_5_ECC_ADJ -3 // largest prob occurs when t2 == t1+this
Prob 
prior5(char t1, char t2, int topLevel) {
    if ((t2 <= t1 + PRIOR_5_ECC_ADJ - PRIOR_4_L) || (t2 >= t1 + PRIOR_5_ECC_ADJ + PRIOR_4_M)) 
        return (Prob)PRIOR_5_FLOOR;
    else {
        if (t2 <= t1+PRIOR_5_ECC_ADJ) {
            return (Prob)topLevel - (Prob)(topLevel - PRIOR_5_FLOOR) * (Prob) (t1 + PRIOR_5_ECC_ADJ - t2) / (Prob)PRIOR_5_L ;
        } else {
            return (Prob)topLevel - (Prob)(topLevel - PRIOR_5_FLOOR) * (Prob) (t2 - t1 - PRIOR_5_ECC_ADJ) / (Prob)PRIOR_5_M ;
        }
    }
}//prior5()

/*
** Read weights, one per line, into prior.
** Normalisation not performed.
** Checkes that there is the right number of weights.
*/
void
readPriorFile(const char *filename, Prob *prior) {
   #define BUFF_SIZE 128
   char buff[BUFF_SIZE];
   FILE *file = fopen(filename, "r");
   if (file == NULL) {
      fprintf(stderr,"Cannot open file %s for reading\n",filename);
      exit(-1);
   }

   int i = 0;
   while (fgets(buff, BUFF_SIZE, file) != NULL) {
      if (i >= domain) {
         fprintf(stderr,"Too many lines in %s. Expected %d.\n",filename,domain);
         exit(-1);
      }
      sscanf(buff, "%f", &prior[i++]);
   }

   fclose(file);
}//readPriorFile()


/*
** Create the prior so it sits around 
** p == 1 is uniform prior
** p == 2 is "tri level" (same/diff-by-1/diff-by-many)
*/
void
makePrior(int p) {
    prior = (Prob *)smalloc(sizeof(Prob) * domain);
    if (p == 0) {
        printf("\n*********************************************************************\n");
        printf("WARNING: Be careful that your tree file is correct when using prior==0!");
        printf("\n*********************************************************************\n");
        readPriorFile(PRIOR_FILENAME, prior);
    } else if (p == 1) {
        for(int i = 0 ; i < domain ; i++)
            prior[i] = (Prob)1 / (Prob)domain;
    } else if (p == 2) {
        for(int i = 0 ; i < domain ; i++) {
            char *ts = digits[i];
            char maxDb = ts[0], minDb=ts[0];
            for(int j = 1 ; j < numberOfLocations ; j++) {
                if (ts[j] < minDb)
                    minDb = ts[j];
                if (ts[j] > maxDb)
                    maxDb = ts[j];
            }

               if (maxDb - minDb > 1)
                  prior[i] = (Prob)1;
               else if (maxDb - minDb == 1)
                  prior[i] = (Prob)2;
               else 
                  prior[i] = (Prob)4;
        }
    } else if (p == 3) { 
        // HAS some zero values!
        // p==3 == skew normal m=t1, s=5, alpha=-4 for each subdomain
        // only works for 2 locations
        // f <- function(t2,a,t1,s) { 2/s* dnorm((t2-t1)/s)*pnorm(a*(t2-t1)/s) }
        if (numberOfLocations != 2) {
            fprintf(stderr, "Trying to use prior==3 for n!=2, aborting\n");
            exit(-1);
        }
        #define PRIOR_3_SIGMA 10.0
        #define PRIOR_3_ALPHA -4.0
        for(int i = 0 ; i < domain ; i++) {
            char t1 = digits[i][0];
            char t2 = digits[i][1];
            float t2_standadised = ((float)t2 - (float)t1) / PRIOR_3_SIGMA;
            double normPdf = gsl_ran_ugaussian_pdf(t2_standadised);
            double normCdf = gsl_cdf_ugaussian_Q  (PRIOR_3_ALPHA * t2_standadised);
            prior[i] = (Prob)(2.0/PRIOR_3_SIGMA * normPdf * normCdf);
        }
    } else if (p == 4) { 
        if (numberOfLocations != 2) {
            fprintf(stderr, "Trying to use prior==4 for n!=2, aborting\n");
            exit(-1);
        }
        for(int i = 0 ; i < domain ; i++) {
            char t1 = digits[i][0];
            char t2 = digits[i][1];
            prior[i] = prior4(t1, t2);
        }
    } else if ((p >= 5) && (p <= 7)) { 
        if (numberOfLocations != 2) {
            fprintf(stderr, "Trying to use prior==5 for n!=2, aborting\n");
            exit(-1);
        }
        for(int i = 0 ; i < domain ; i++) {
            char t1 = digits[i][0];
            char t2 = digits[i][1];
            prior[i] = prior5(t1, t2, (p == 5) ? 1 : 50) * (p == 7 ? gsl_ran_ugaussian_pdf(((float)t1 - 31.0)/20.0) : 1.0);
        }
    } else if ((p >= 20) && (p <= 60)) {
        if (numberOfLocations != 1) {
            fprintf(stderr, "Trying to use biModal pdf prior for n!=1, aborting\n");
            exit(-1);
        }
        float gpdf[] = { 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.2, 0.3, 0.2, 0.15,0.1, 0.02,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001 };
      float npdf[] = { 0.001, 0.009, 0.03, 0.05, 0.1, 0.2, 0.3, 0.2, 0.05, 0.025, 0.01, 0.001};
        float newPdf[domain];

        for(int i = 0 ; i < domain ; i++)
            newPdf[i] = 0.001;
        
        int mode = p-20;
        for(int i = 0 ; i < 12 ; i++)
            if ((mode-6+i >= 0) && (mode-6+i < domain))
                newPdf[mode-6+i] = npdf[i];

        for(int i = 0 ; i < domain ; i++) 
            prior[i] = gpdf[i] + 4*newPdf[i];
    }

        // normalise for all but p==1
    if (p != 1) {
        Prob sum = (Prob)0;
        for(int i = 0 ; i < domain ; i++)
            sum += prior[i];
        for(int i = 0 ; i < domain ; i++)
{
            prior[i] /= sum;
//if (prior[i] < 0) {
//char t1 = digits[i][0];
//char t2 = digits[i][1];
//char t3 = digits[i][2];
//printf("%2d %2d %2d %20.18f\n",t1,t2,t3,prior[i]);
//}
}
    }
}//makePrior()


/*
** Set global variables from command line.
** Load tree if it exists.
** RETURN 0 on success
*/
char
processCommandLine(int argc,char *argv[]) {
    sscanf(argv[1], "%d", &numberOfTrials);
    sscanf(argv[2], "%d", &numberOfLocations);
    sscanf(argv[3], "%d", &numberOfDbValues);

    int i;
    dbValues = (int *)smalloc(sizeof(int)*numberOfDbValues);
    for(i = 4 ; i < 4 + numberOfDbValues; i++)
        sscanf(argv[i], "%d", dbValues + i - 4);

    fos = (Fos *)smalloc(sizeof(Fos)*numberOfLocations);
    int j;
    int p = 0;
    for( j = i ; j < i + 4*numberOfLocations ; p++) {
        sscanf(argv[j], "%f", &(fos[p].threshold)); j++;
        sscanf(argv[j], "%f", &(fos[p].slope)); j++;
        sscanf(argv[j], "%f", &(fos[p].fp)); j++;
        sscanf(argv[j], "%f", &(fos[p].fn)); j++;
    }

    sscanf(argv[j]  , "%d", &priorType);
    sscanf(argv[j+1], "%d", &lfType);
    sscanf(argv[j+2], "%c", &stopType);
    sscanf(argv[j+3], "%f", &stopValue);
    sscanf(argv[j+4], "%d", &maxDepth);

    rewriteTreeFlag = argv[j+5][0] == 'w';

    verbose = 0;
    if (j+6 < argc)
        verbose = argv[j+6][0] == 'v';

        // now validate a few things

    if (stopType == 'H')
        stopFunction = stopEntropy;
    else if (stopType == 'N')
        stopFunction = stopNumPres;
    else if (stopType == 'S')
        stopFunction = stopStdev;
    else {
        fprintf(stderr,"Invalid stopType %c\n", stopType);
        return 1;
    }

    if ((priorType > 7) && ((priorType <20) || (priorType > 60))) {
        fprintf(stderr,"Unknown prior type %d\n",priorType);
        return 1;
    }
    if (lfType != 1) {
        fprintf(stderr,"Unknown lf type %d\n",lfType);
        return 1;
    }

    domain = (int)pow(numberOfDbValues, numberOfLocations);

        // these are the digits of the base |D| expansion of a domain value
    digits = (char **)smalloc(sizeof(char *) * domain);
    digits[0] = (char *)smalloc(sizeof(char) * numberOfLocations);
    for(int j = 0 ; j < numberOfLocations ; j++)
        digits[0][j] = 0;
    for(int i = 1 ; i < domain ; i++) {
        digits[i] = (char *)smalloc(sizeof(char) * numberOfLocations);
        for(int j = 0 ; j < numberOfLocations ; j++)
            digits[i][j] = digits[i-1][j];

            // increment the rightmost digit, which might
            // cascade to other digits (watch out for digit[-1])
        int d = numberOfLocations-1;
        digits[i][d]++;
        while (d >= 0 && digits[i][d] >= numberOfDbValues) {
            digits[i][d] = 0;
            if (d > 0)
                digits[i][d-1]++;
            d--;
        }
    }

    makePrior(priorType);
//for(int i = 0 ; i < domain ; i++)
//printf("%10d %f\n",i,prior[i]);
    makeLf(lfType);

        // Check if there is a tree file for this config
        // If so, load it, else create an empty tree and 
        // put prior into slot "0" ready for first use.
    tree = treeMakeOrLoad(numberOfLocations, numberOfDbValues, dbValues, priorType, lfType, stopType, stopValue, maxDepth);

    return 0;
}//processCommandLine()

/*
** Set values in m (mean, sd, entropy) of pdf.
*/
void
moments(Prob *pdf, MomentType *m) {
    Prob sum=0, sumSq=0, H=0;
    for(int i = 0 ; i < domain ; i++) 
        if (pdf[i] > 0) {
            sum   += pdf[i] * i;
            sumSq += pdf[i] * i * i;
            H     += ENTROPY(pdf[i]);
        }

    m->mean = sum;
    m->sd   = sqrt(sumSq - sum*sum);
    m->H    = H;
}// moments()

/*
** Work back up the tree from tp to root making a stack of node ptrs, 
** then work back down again making the pdf from the prior.
** Stores pdf in pdf, 
** ASSUMES it is non-NULL.
*/
void
makePdf(uint tp, Prob *pdf) {
//fprintf(stderr,"Make pdf for node %d\n",tp);
    TreeNode *path[1024];        // ASSUMES no path will be longer than this
    int nextP = 0;
    while (tp != TREE_NULL) {
        TreeNode *p = TO_N_PTR(tp);
        path[nextP++] = p;
//fprintf(stderr,"Make pdf : added %d to the stack %lx\n",tp,p);
        tp = p->parent;
    } 
    path[nextP] = TO_N_PTR(0);    // because TREE_NULL is root is 0, need to explicitly add root
//fprintf(stderr,"Make pdf : added %d to the stack %lx\n",0,path[nextP]);

    for(int i = 0 ; i < domain ; i++)
        pdf[i] = prior[i];

//fprintf(stderr,"Make pdf  == prior\n");
    for( ; nextP > 0 ; nextP--) {
        TreeNode *p = path[nextP];
        TreeNode *child = path[nextP-1];

        int location = p->location, dbIndex = p->dbIndex;
        Prob sum = 0.0;
        if (child->type == TREE_NODE_TYPE_NO_PDF) {
//fprintf(stderr,"Make pdf multiply %lx * NOLF %d %d \n", p, location, dbIndex);
            for(int i = 0 ; i < domain ; i++) {
                pdf[i] *= (1.0 - lfLut[dbIndex][ (int)digits[i][location] ]);
                sum += pdf[i];
            }
        } else {
//fprintf(stderr,"Make pdf multiply %lx * YES LF %d %d \n", p, location, dbIndex);
            for(int i = 0 ; i < domain ; i++) {
                pdf[i] *= lfLut[dbIndex][ (int)digits[i][location] ];
                sum += pdf[i];
            }
        }

        for(int i = 0 ; i < domain ; i++) 
            pdf[i] /= sum;
    }
}//makePdf()

/*
** Return true if stop, false otherwise
*/
char stopStdev(Prob *pdf, int depth)   { MomentType m; moments(pdf, &m); return m.sd <= stopValue; }
char stopEntropy(Prob *pdf, int depth) { MomentType m; moments(pdf, &m); return m.H <= stopValue;  }
char stopNumPres(Prob *pdf, int depth) { return depth >= stopValue; }

/*
** Compute the expected set of db values, one for each location.
** ASSUMES answer has room [0..numLocations-1].
*/
void
expectedSample(Prob *pdf, float *answer) {
    for(int i = 0 ; i < numberOfLocations ; i++)
        answer[i] = 0;

    for(int i = 0 ; i < domain ; i++) {
        for(int j = 0 ; j < numberOfLocations ; j++) {
            answer[j] += pdf[i] * dbValues[(int)digits[i][j]];
        }
    }
    return;
}//expectedSample()

/*
** Choose a stimulus
** RETURNS location (0..numberOfLocations-1) and dbIndex (0..numberOfDbValues-1) by reference
*/
void
chooseStimulus(Prob *pdf, int prevLocation, int *location, int *dbIndex) {
    Prob minH = 999999999;
    int minL=-1, minD=-1;
    for(int loc = 0 ; loc < numberOfLocations ; loc++) {
        //if ((numberOfLocations > 1) && (loc == prevLocation))
        //    continue;
        for(int db = 0 ; db < numberOfDbValues ; db++) {

            Prob *yesPdf = workArea1;
            Prob p = 0.0;
            for(int i = 0 ; i < domain ; i++) {
                yesPdf[i] = pdf[i] * lfLut[db][ (int)digits[i][loc] ];
                p += yesPdf[i];
            }
            for(int i = 0 ; i < domain ; i++) 
                yesPdf[i] /= p;

            MomentType m;
            moments(yesPdf, &m);
            Prob EH = p * m.H;

            if (EH < minH) {                // don't need to do more if already too high
                Prob *noPdf = workArea1;     // be careful here as we are overwriting yesPdf via workArea1
                for(int i = 0 ; i < domain ; i++)
                    noPdf[i] = pdf[i] - yesPdf[i] * p;  // no = pdf * (1.0 - lf);
                for(int i = 0 ; i < domain ; i++) 
                    noPdf[i] /= (1-p);
                moments(noPdf, &m);
                EH += (1-p) * m.H;

                if (EH < minH) {
                    minH = EH;
                    minL = loc;
                    minD = db;
                }
            }
        }
    }

    *location = minL;
    *dbIndex  = minD;
}//chooseStimulus()

/*
** Make a new node of either yes or no type at depth depth in the tree.
**    location is location at which stim was presented.
**    dbIndex is index into dbValues at which stim was presented.
**    depth is depth of this node, not parent.
**    pdf is posteria pdf of parent, prior of current.
**    parentPtr is index into nodes chunk of parent
** RETURNS ptr to new node (as uint index into chunk tree->nodes)
*/
uint 
makeNewNode(int location, int dbIndex, Prob *pdf, char nodeType, int depth, uint parentPtr, char printMoments) {
    Prob *newPdf = workArea2;
    Prob p = 0.0;
    if (nodeType == TREE_NODE_TYPE_NO_PDF) {
        for(int i = 0 ; i < domain ; i++) {
            newPdf[i] = pdf[i] * (1.0 - lfLut[dbIndex][ (int)digits[i][location] ]);
            p += newPdf[i];
        }
    } else {
        for(int i = 0 ; i < domain ; i++) {
            newPdf[i] = pdf[i] * lfLut[dbIndex][ (int)digits[i][location] ];
            p += newPdf[i];
        }
    }
    for(int i = 0 ; i < domain ; i++) 
        newPdf[i] /= p;

    if (printMoments) {
        MomentType m;
        moments(newPdf, &m);
        printf("# Moments:  %9.3f %9.3f %9.3f\n",m.mean,m.sd,m.H);
    }
    
    uint ptr;
    TreeNode *node = (TreeNode *)chunkGetNext(tree->nodes, &ptr);
//fprintf(stderr,"Adding new node at %d...",ptr);
//MomentType m;
//moments(newPdf, &m);
//fprintf(stderr,"Moments:  %9.3f %9.3f %9.3f ",m.mean,m.sd,m.H);
    node->valid  = TREE_NODE_VALID;
    node->parent = parentPtr;
    if (depth == maxDepth || stopFunction(newPdf, depth)) {
        uint dataPtr;
        float *temp = (float *)chunkGetNext(tree->leaves, &dataPtr);
        expectedSample(pdf, temp);
        node->type  = TREE_NODE_TYPE_LEAF;
        node->data = dataPtr;
//fprintf(stderr,"Leaf\n");
    } else {
        chooseStimulus(newPdf, location, &location, &dbIndex);

        node->location = location;
        node->dbIndex  = dbIndex;
        node->type     = nodeType;
//fprintf(stderr,"Stim: %d %d type=%s\n",location,dbIndex,nodeType == TREE_NODE_TYPE_NO_PDF ? "no" : "yes");
    }

    return(ptr);
}//makeNewNode()

/*
** Loop until stop met or maxDepth, returning (by reference)
** number of presentaions per location and (by value) expected value per location.
** ASSUMES n has already been malloced.
** ASSUMES tree has a valid root
**
** Need to be really careful using chunkGetNext as it can realloc, making old
** pointers into chunks invalid. That's why TO_N_PTR is sprinkled everywhere currently.
*/
float *
buss(int *n) {
    for(int j = 0 ; j < numberOfLocations ; j++)    // num presentations per location
        n[j] = 0;

    int depth = 0;
    uint currTreeNode = treeGetRoot(tree);
    for(;;) {
        // assert: now currTreeNode is a valid node
        TreeNode *currTptr = TO_N_PTR(currTreeNode);
//printf("currTptr = TO_(%d)\n",currTreeNode);fflush(stdout);
        if (currTptr->type == TREE_NODE_TYPE_LEAF) {
            return (float *)chunkGet(tree->leaves, currTptr->data);
        } else {
            int location = currTptr->location, 
                dbIndex  = currTptr->dbIndex;
            
            if (verbose)
                printf("# Stim %03d %2d %2d ",depth, location, dbValues[dbIndex]);
            
            //char seen = (double)random()/(double)((1 << 31)-1) <= PSI(dbValues[dbIndex], fos[location]);
            char seen = gsl_rng_uniform(rng) <= PSI(dbValues[dbIndex], fos[location]);

            n[location]++;
            depth++;

            if (verbose)
                printf("%s\n", seen ? "Yes" : "No");
            
                // add both valid yes and no nodes
            if ((currTptr->no == TREE_NULL) || (currTptr->yes == TREE_NULL)) {
                Prob *pdf = workArea3;
                makePdf(currTreeNode, pdf);
            
                if (currTptr->no == TREE_NULL) {      // note makeNewNode might move currTreeNode in mem, so use temp = ...
                    uint temp = makeNewNode(location, dbIndex, pdf, TREE_NODE_TYPE_NO_PDF, depth, currTreeNode, verbose &&
!seen);
                    TO_N_PTR(currTreeNode)->no = temp;
                }

                if (TO_N_PTR(currTreeNode)->yes == TREE_NULL) {
                    uint temp = makeNewNode(location, dbIndex, pdf, TREE_NODE_TYPE_YES_PDF, depth, currTreeNode, verbose &&
seen);
                    TO_N_PTR(currTreeNode)->yes = temp;
                }

                treeChanged = 1;
            }

            currTreeNode = seen ? TO_N_PTR(currTreeNode)->yes : TO_N_PTR(currTreeNode)->no;
        }
    }
}// buss()

int
main(int argc, char *argv[])  {

    if ((argc < 15) || processCommandLine(argc, argv)) {
        fprintf(stderr, "Usage %s num-trials num-locs num-dB {dBs}+ {true-thresh slope fp fn}+ priorFlag lfFlag stopType stopVal maxDepth {w|n} [v]\n",argv[0]);
        fprintf(stderr, "   where\n");
        fprintf(stderr, "       priorNum = 0 - Read prior from file %s\n", PRIOR_FILENAME);
        fprintf(stderr, "                = 1 - uniform pdf\n");
        fprintf(stderr, "                = 2 - tri-level pdf with abs diff in thresh {>1, 1, 0} weighted as {1,2,4}\n");
        fprintf(stderr, "                = 3 - skew normal pdf - warning has some probs == 0 \n");
        fprintf(stderr, "                = 4 - quad-linear with dynamic peak (hence (0,0) dominates)!\n");
        fprintf(stderr, "                = 5 - tri-linear with fixed peak of 2!\n");
        fprintf(stderr, "                = 6 - tri-linear with fixed peak of 50!\n");
        fprintf(stderr, "                = 7 - as for 6 with Loc1 biased with Gauss(31,20)\n");
        fprintf(stderr, "                = 20..60 - bimodal pdf with normal mode at value-20\n");
        fprintf(stderr, "                     (eg 20==0dB, 30==10dB)\n");
        fprintf(stderr, "       lfNum    = 1 - Henson -0.098 & 3.62 capped at 1.0, fp=fn=3%%\n");
        fprintf(stderr, "    stopType    = H | N | S\n");
        fprintf(stderr, "       {w|n} w == rewrite tree file at end if tree changed, n == do not rewrite tree file\n");
        fprintf(stderr, "       v == verbose (optional) \n");
        return -1;
    }

    rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, time(NULL));

    //srandom(time(NULL));
    printf("\n# Trials     = %d",numberOfTrials);
    printf("\n# |L|        = %d",numberOfLocations);
    printf("\n# D (size=%2d)= {",numberOfDbValues);
    for(int i = 0 ; i < numberOfDbValues; i++)
        printf(" %d", dbValues[i]);
    printf("}");
    for(int i = 0 ; i < numberOfLocations; i++)
        printf("\n# FOSS %2d    = %4.2f %4.2f %4.2f %4.2f",i,fos[i].threshold, fos[i].slope, fos[i].fp, fos[i].fn);
    printf("\n# Prior      = %d",priorType);
    printf("\n# Lf         = %d",lfType);
    printf("\n# Stop type  = %c",stopType);
    printf("\n# Stop val   = %f",stopValue);
    printf("\n# Max Depth  = %d",maxDepth);
    printf("\n# Domain     = %d",domain);
    printf("\n");

        // for chooseStimulus
    workArea1 = (Prob *)smalloc(sizeof(Prob) * domain);
    workArea2 = (Prob *)smalloc(sizeof(Prob) * domain);
    workArea3 = (Prob *)smalloc(sizeof(Prob) * domain);

        // if the tree doesn't have a valid root, make one
    if (!treeHasValidRoot(tree)) { 
        Prob *pdf = prior;
        int location, dbIndex;
        if (verbose) {
            MomentType m;
            moments(pdf, &m);
            printf("# Moments:  %9.3f %9.3f %9.3f\n",m.mean,m.sd,m.H);
        }
        chooseStimulus(pdf, -1, &location, &dbIndex);
        treeCreateRoot(tree, location, dbIndex);
        treeChanged = 1;
    } else 
        treeChanged = 0;

        // now call buss per trial and output
    int   *numPres = (int *)smalloc(sizeof(int) * numberOfLocations);
    for(int i = 0 ; i < numberOfTrials ; i++) {
        struct timeval tp1, tp2;

        gettimeofday(&tp1, NULL);

        float *answer = buss(numPres);

        gettimeofday(&tp2, NULL);

        printf("%4d ",i);
        for(int j = 0 ; j < numberOfLocations ; j++)
            printf(" %2d",numPres[j]);
        printf("   ");
        for(int j = 0 ; j < numberOfLocations ; j++)
            printf(" %5.2f",answer[j]);

        double s1 = (double)tp1.tv_sec + (double)tp1.tv_usec/1000000.0;
        double s2 = (double)tp2.tv_sec + (double)tp2.tv_usec/1000000.0;
        printf("   %5.1f seconds\n", s2 - s1);

    }

    free(dbValues);
    free(fos);
    free(numPres);
    free(workArea1);
    free(workArea2);
    free(workArea3);
    gsl_rng_free (rng);
    //free(digits);  // bad! memory leak
    //free(lfLut);   // bad! memory leak

    if (treeChanged && rewriteTreeFlag) {
        fprintf(stderr,"# Writing tree to file...\n");
        treeSave(tree); // MEMORY LEAK - should free it
    }

    //report_memuse();

    return 0;
}
