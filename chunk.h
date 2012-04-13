/*
** Routines to handle one big array that is nominally an 
** array of elements of size "size".
**
** Andrew Turpin
** Wed  9 Mar 2011 13:19:08 EST
*/
typedef unsigned int uint;
typedef unsigned char uchar;

typedef struct chunk {
    uint size;  // number of bytes in each "element" of chunk
    uint max;   // p[0..(max-1)*size] is a valid index
    uint n;     // p[0..(n-1)*size] is used
    uchar *p;   // one big array 
} Chunk;

#define MAX_ELEMENTS_TO_REALLOC (1 << 10)
#define min(_a,_b) ((_a) < (_b) ? (_a) : (_b))

Chunk *chunkCreate(uint size, uint max, void (*init)(void *));
void *chunkGetNext(Chunk *c, uint *index);   // get address of p[n*size] and n in *index, increment c->n
uint chunkGetN(Chunk *c);          // get number of elements in chunk
void *chunkGet(Chunk *c, uint i);  // get ptr to element i in chunk
void chunkWrite(Chunk *c, int f);
Chunk *chunkRead(int fd);
