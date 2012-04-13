// smalloc.c
// a safer, though simple, version of malloc and realloc
// Anthony Wirth
// June 2010

#include <stdlib.h>
#include <stdio.h>
#include "smalloc.h"

static size_t m_counter = 0;
static size_t total_m = 0;
static size_t max_m = 0;
static size_t w_counter = 0;
static size_t r_counter = 0;

void *
smalloc(size_t size)
{
	void *p = malloc(size);
	if(p == NULL){
		fprintf(stderr,"There was not enough space on the heap. Aborting.\n");
		exit(EXIT_FAILURE);
	}
	m_counter += size;
	if(m_counter > max_m){
		max_m = m_counter;
	}
	total_m += size;
	return p;
}

void *
srealloc(void *p, size_t size, size_t old_size)
{
	void *q;
	q = realloc(p,size);
	if(q == NULL){
		fprintf(stderr,"There was not enough space on the heap. "
		"Aborting.\n");
		exit(EXIT_FAILURE);
	}
	m_counter += (size-old_size);
	if(m_counter > max_m){
		max_m = m_counter;
	}
	total_m += (size-old_size);
	return q;
}

void
sfwrite(const void *ptr, size_t size, size_t nitems,
	FILE *stream)
{
	size_t written;
	w_counter += size * nitems;
	written = fwrite(ptr,size,nitems,stream);
	if(written < nitems){
		fprintf(stderr,"Was not able to write all bytes to file. "
			"Aborting.\n");
		exit(EXIT_FAILURE);
	}
}

size_t
sfread(void *ptr, size_t size, size_t nitems,  FILE *stream)
{
	size_t nread;
	nread = fread(ptr,size,nitems,stream);
	r_counter += nread*size;
	return nread;
}

void
sfree(void *p,size_t size)
{
	free(p);
	m_counter -= size;
}

void
report_memuse(void)
{
	fprintf(stderr,"# %12ld bytes total allocated.\n",total_m);
	fprintf(stderr,"# %12ld bytes max at once allocated.\n",max_m);
	fprintf(stderr,"# %12ld bytes written.\n",w_counter);
	fprintf(stderr,"# %12ld bytes read.\n",r_counter);
}
