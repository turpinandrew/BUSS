#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "smalloc.h"
#include "chunk.h"

Chunk *
chunkCreate(uint size, uint max, void (*init)(void *)) {
    Chunk *c = (Chunk *)smalloc(sizeof(Chunk));
    c->p  = (uchar *)smalloc(size * max);
    c->n = 0;
    c->max = max;
    c->size = size;

        // call the init function on each element if it exists
    if (init != NULL) {
        for(uint i = 0 ; i < max ; i++)
            init( c->p + c->size * i);
    }
    
    return c;
}//chunkCreate()

/*
** Return the address of p[n*size] element and the value of n by ref in *i.
** Increment n. 
** if n >= max realloc first so p[n*size] is valid.
*/
void *
chunkGetNext(Chunk *c, uint *i) {
    if (c->n >= c->max) {
        int extraElements = min(c->max, MAX_ELEMENTS_TO_REALLOC);
        c->p = (uchar *)srealloc(c->p, (c->max + extraElements) * c->size, c->max * c->size);
        c->max += extraElements;
    }
    *i = c->n;
    c->n += 1;
    return chunkGet(c, *i);
}//safeGet()

uint chunkGetN(Chunk *c) { return c->n; }

void *chunkGet(Chunk *c, uint i) { return c->p + i * c->size; } // no valid check

void 
chunkWrite(Chunk *c, int f) {
    write(f, c, sizeof(Chunk));
    write(f, c->p, sizeof(uchar) * c->n * c->size);
}//chunkWrite()

Chunk *
chunkRead(int fd) {
    Chunk *c = (Chunk *)smalloc(sizeof(Chunk));
    read(fd, c,    sizeof(Chunk));
    c->p = (uchar *)smalloc(c->size  * c->n);
    read(fd, c->p, sizeof(uchar) * c->n * c->size);
    c->max = c->n;

    return c;
}//chunkRead()
