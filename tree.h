/*
** Decision tree for loading/saving, etc
** Pointers all done as indexes into chunks to allow
** realloc of chunks, and easy serialisation.
** 
** Wed  9 Mar 2011 10:31:25 EST
*/

#define MAX_PRIOR_FILENAME_LENGTH 512

typedef unsigned char uchar;
typedef unsigned int  uint;

typedef struct treeNode { 
    uchar location;
    uchar dbIndex;
    uchar type;
    uchar valid;
    uint data;     // could be a pdf (index into pdfs chunk) or a leaf (index into leaves chunk)
    uint yes, no, parent;  // ptrs as indexes into nodes[] chunk
} TreeNode;

    // Tree consists of two blocks of memory
    // nodes[], leaves[] 
    // The tree structure is threaded through nodes[] with data pointers
    // into leaves[].
    // This makes load/save easier.
typedef struct tree {
    uint numLocations;  // stored so we know length of leaves
    uint numDbValues;   // stored so we can calc length of pdfs

    Chunk *nodes;
    Chunk *leaves;

    char priorType;     // so we can validate on load
    char lfType;        
    char stopType;
    double stopValue;
    int maxDepth;
    int *dbValues;
} Tree;

#define TREE_INITIAL_NUM_NODES  (1 << 9)
#define TREE_INITIAL_NUM_LEAVES (1 << 9)
#define TREE_INITIAL_NUM_PDFS   (1 << 9)

#define TREE_NULL ((uint)0)  // index of NULL in nodes[]. 0 is the index of the root, but we are reusing it as NULL

#define TREE_NODE_INVALID 0
#define TREE_NODE_VALID   1

#define TREE_NODE_TYPE_ROOT    0
#define TREE_NODE_TYPE_NO_PDF  1
#define TREE_NODE_TYPE_YES_PDF 2
#define TREE_NODE_TYPE_LEAF    4

Tree *treeMakeOrLoad(uint n, uint d, int *dbValues, char pType, char lfType, char stopT, double stopV, int maxD);
int treeHasValidRoot(Tree *t);
void treeCreateRoot(Tree *t, int location, int dbIndex);
uint treeGetRoot(Tree *t);
void treeSave(Tree *t);
Tree *treeLoad(char *filename);
