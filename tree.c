/*
** Tree of decision nodes.
** Root is in tree->nodes[0].
**
** Andrew Turpin
** Wed  9 Mar 2011 11:28:57 EST
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>
#include "smalloc.h"
#include "chunk.h"
#include "tree.h"

// helper for makeTree()
void initNode(void *pp) { 
   TreeNode *p = (TreeNode *)pp; 
   p->valid = TREE_NODE_INVALID; 
   p->yes = p->no = p->parent = TREE_NULL; 
}

/*
** Allocate all memory, etc
*/
Tree *
makeTree(int n, int d, int *dbValues, char pType, char lfType, char stopT, double stopV, int maxD) {
    Tree *t   = (Tree *)smalloc(sizeof(Tree));

    t->numLocations = n;
    t->numDbValues  = d;

    t->nodes  = chunkCreate(sizeof(TreeNode),  TREE_INITIAL_NUM_NODES , initNode);
    t->leaves = chunkCreate(sizeof(float) * n, TREE_INITIAL_NUM_LEAVES, NULL);

    t->priorType  = pType;     
    t->lfType     = lfType;        
    t->stopType   = stopT;
    t->stopValue  = stopV;
    t->maxDepth   = maxD;
    t->dbValues = (int *)smalloc(sizeof(int) * d);
    for(int i = 0 ; i < d ; i++)
        t->dbValues[i] = dbValues[i];

    return t;
}//makeTree()

/*
** Create a filename for the tree that is hopefully
** unique to this procedure.
*/
char *
makeFilename(int n, int d, int *dbValues, char pType, char lfType, char stopT, double stopV, int maxD) {
    char buff[32];
    char *filename = (char *)smalloc(sizeof(char)*4 * 1024);  // just to be safe
    filename[0] = 0;
    
    strcat(filename, "tree");
    sprintf(buff, "_%d", n);      strcat(filename, buff);
    sprintf(buff, "_%d", d);      strcat(filename, buff);
    sprintf(buff, "_%d", pType);  strcat(filename, buff);
    sprintf(buff, "_%d", lfType); strcat(filename, buff);
    sprintf(buff, "_%c", stopT);  strcat(filename, buff);
    sprintf(buff, "_%f", stopV);  strcat(filename, buff);
    sprintf(buff, "_%d", maxD);   strcat(filename, buff);

    for(int i = 0 ; i < d ; i++) {
        sprintf(buff, "_%d", dbValues[i]);   
        strcat(filename, buff);
    }

    strcat(filename, ".bin");

    return filename;
}//makeFilename()

/*
** Attempt to load the tree, otherwise initialise an empty tree
*/
Tree *
treeMakeOrLoad(uint n, uint d, int *dbValues, char pType, char lfType, char stopT, double stopV, int maxD) {
    Tree *t;
    char *filename = makeFilename(n, d, dbValues, pType, lfType, stopT, stopV, maxD);
    FILE *f = fopen(filename, "r");
    if (f == NULL)
        t = makeTree(n, d, dbValues, pType, lfType, stopT, stopV, maxD);
    else {
        fclose(f);
        t = treeLoad(filename);

        #define CHECK(_t, _a,_b) if ((_a) != (_b)) fprintf(stderr,"ERROR reading tree: %s %d != %d\n", _t, _a, _b);
        #define CHECKF(_t, _a,_b) if ((_a) != (_b)) fprintf(stderr,"ERROR reading tree: %s %f != %f\n", _t, _a, _b);
        CHECK("numL",t->numLocations,n) 
        CHECK("numD",t->numDbValues,d) 
        CHECK("pType",t->priorType,pType) 
        CHECK("maxD",t->maxDepth,maxD) 
        CHECK("stopT",t->stopType,stopT) 
        CHECKF("stopV",t->stopValue,stopV) 
        for(uint i = 0 ; i < d ; i++)
            CHECK("dbV",t->dbValues[i], dbValues[i]);
    }

    return t;
}//treeMakeOrLoad()

/*
** Return pointer to root of tree (nodes[0])
*/
uint
treeGetRoot(Tree *t) {
    if (t == NULL || chunkGetN(t->nodes) < 1)
        fprintf(stderr,"ERROR: trying to get root of tree before calling treeMakeOrLoad(...)\n");
    return 0;
}//treeGetRoot()

/*
** return 1 if tree root is valid, 0 otherwise
*/
int
treeHasValidRoot(Tree *t) {
    if (t == NULL || chunkGetN(t->nodes) < 1)
        return 0;
    TreeNode *root = (TreeNode *)chunkGet(t->nodes, treeGetRoot(t));
    return root->type == TREE_NODE_TYPE_ROOT && root->valid != TREE_NODE_INVALID;
}//treeHasValidRoot()

/*
** Initialise tree->nodes[0].
** and set tree->numTreeNodes == 1
*/
void 
treeCreateRoot(Tree *t, int location, int dbIndex) {
    if (t == NULL) {
        fprintf(stderr,"ERROR: called treeCreateRoot(...) before calling treeMakeOrLoad(...).\n");
        return;
    }

    uint index;
    TreeNode *n = (TreeNode *)chunkGetNext(t->nodes, &index);

    if (index != 0) {
        fprintf(stderr,"ERROR: called treeCreateRoot(...) without an empty tree. The root must be in position 0\n");
        fprintf(stderr,"       of the nodes chunk or disaster will ensue!\n");
        return;
    }

    n->location = (uchar)location;
    n->dbIndex  = (uchar)dbIndex;
    n->valid    = TREE_NODE_VALID;
    n->type     = TREE_NODE_TYPE_ROOT;
    //n->data     = NULL;                 // not required as NULL is meaningless here
    n->yes      = TREE_NULL;
    n->no       = TREE_NULL;
    n->parent   = TREE_NULL;
}//treeCreateRoot()

void
treeSave(Tree *t) {
    char *filename = makeFilename(t->numLocations, t->numDbValues, t->dbValues, t->priorType, 
                                  t->lfType, t->stopType, t->stopValue, t->maxDepth);
    //struct flock fl;
    //fl.l_type   = F_WRLCK;  /* F_RDLCK, F_WRLCK, F_UNLCK    */
    //fl.l_whence = SEEK_SET; /* SEEK_SET, SEEK_CUR, SEEK_END */
    //fl.l_start  = 0;        /* Offset from l_whence         */
    //fl.l_len    = 0;        /* length, 0 = to EOF           */
    //fl.l_pid    = getpid(); /* our PID                      */
    
    int f = open(filename, O_WRONLY | O_CREAT); // , S_IRWXU);
    if (f == -1) {
        perror("open tree file:");
        return;
    }
    //fcntl(f, F_SETLKW, &fl);  /* set the lock, waiting if necessary */
    write(f, t, sizeof(Tree));
    write(f, t->dbValues, sizeof(int) * t->numDbValues);
    chunkWrite(t->nodes,    f);
    chunkWrite(t->leaves,   f);
    //fl.l_type   = F_UNLCK;  /* tell it to unlock the region */
    //fcntl(f, F_SETLK, &fl); /* set the region to unlocked */
    close(f);
}//saveTree()

Tree *
treeLoad(char *filename) {
    //struct flock fl;
    //fl.l_type   = F_RDLCK;  /* F_RDLCK, F_WRLCK, F_UNLCK    */
    //fl.l_whence = SEEK_SET; /* SEEK_SET, SEEK_CUR, SEEK_END */
    //fl.l_start  = 0;        /* Offset from l_whence         */
    //fl.l_len    = 0;        /* length, 0 = to EOF           */
    //fl.l_pid    = getpid(); /* our PID                      */

    int fd = open(filename, O_RDONLY);  /* get the file descriptor */
    //fcntl(fd, F_SETLKW, &fl);  /* set the lock, waiting if necessary */

    Tree *t   = (Tree *)smalloc(sizeof(Tree));
    read(fd, t, sizeof(Tree));
    t->dbValues = (int *)smalloc(sizeof(int) * t->numDbValues);
    read(fd, t->dbValues, sizeof(int) *  t->numDbValues);
    t->nodes  = chunkRead(fd);
    t->leaves = chunkRead(fd);

    //fl.l_type   = F_UNLCK;  /* tell it to unlock the region */
    //fcntl(fd, F_SETLK, &fl); /* set the region to unlocked */
    close(fd);

    return t;
}//treeLoad()
