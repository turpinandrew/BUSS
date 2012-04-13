/*
** Print out the header and a path in a tree file.
**
** Andrew Turpin
** Wed  9 Mar 2011 21:27:44 EST
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "chunk.h"
#include "tree.h"
#include "ent.h"

typedef float Prob; // in case I want double later...

#define MAX_SEQ_LEN 128

    // command line parameters
char *filename;
char sequence[MAX_SEQ_LEN];

int domain;   // |D|^|L| = numberOfDbValues^numberOfLocations

Tree *tree;

int
main(int argc, char *argv[])  {

    if (argc < 3) {
        fprintf(stderr, "Usage %s filename sequence\n", argv[0]);
        fprintf(stderr, "   where\n");
        fprintf(stderr, "       sequence is like YNYYNYNN\n");
        return -1;
    }
    FILE *f = fopen(argv[1], "r");
    if (f == NULL) {
        fprintf(stderr, "Cannot open file %s\n",argv[1]);
        return -1;
    }
    fclose(f);

    tree = treeLoad(argv[1]);
    fclose(f);

    printf("\n# |L|        = %d",tree->numLocations);
    printf("\n# D (size=%2d)= {",tree->numDbValues);
    for(uint i = 0 ; i < tree->numDbValues; i++)
        printf(" %d", tree->dbValues[i]);
    printf("}");
    printf("\n# Prior      = %d",tree->priorType);
    printf("\n# Lf         = %d",tree->lfType);
    printf("\n# Stop type  = %c",tree->stopType);
    printf("\n# Stop val   = %f",tree->stopValue);
    printf("\n# Max Depth  = %d",tree->maxDepth);
    printf("\n");
    printf("\n# Node chunk (size= %d ,max= %d ,n=%d ,p=%lx)",tree->nodes->size,tree->nodes->max,tree->nodes->n,tree->nodes->p);
    printf("\n# Leaf chunk (size= %d ,max= %d ,n=%d ,p=%lx)",tree->leaves->size,tree->leaves->max,tree->leaves->n,tree->leaves->p);

    domain = pow(tree->numDbValues, tree->numLocations);

    char *response = argv[2];
    uint treeNode = treeGetRoot(tree);
    while (*response != 0) {
        printf("Treenode= %10d: ",treeNode);
        TreeNode *t = (TreeNode *)chunkGet(tree->nodes, treeNode);
        if (t->valid == TREE_NODE_INVALID) {
            printf("INVALID\n");
            *response = 0;
        } else {
            if (t->type == TREE_NODE_TYPE_LEAF) {
                float *d = (float *)chunkGet(tree->leaves, t->data);
                printf("Leaves: ");
                for(uint i = 0 ; i < tree->numLocations ; i++)
                    printf("%f ", *(d + i));
                printf("\n");
            } else {
                printf("Stim: %2d %2d Resp: %c\n",t->location, tree->dbValues[t->dbIndex], *response);

                if (*response == 'N')
                    treeNode = t->no;
                else
                    treeNode = t->yes;
            }
            response++;
        }
    }
    TreeNode *t = (TreeNode *)chunkGet(tree->nodes, treeNode);
    if (t->type == TREE_NODE_TYPE_LEAF) {
        float *d = (float *)chunkGet(tree->leaves, t->data);
        printf("Leaves: ");
        for(uint i = 0 ; i < tree->numLocations ; i++)
            printf("%f ", *(d + i));
        printf("\n");
    } 

    return 0;
}
