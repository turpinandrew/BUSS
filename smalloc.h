// smalloc.h
// Anthony Wirth
// June 2010

void
*smalloc(size_t size);

void *
srealloc(void *p, size_t size, size_t old_size);

void
sfwrite(const void *ptr, size_t size, size_t nitems,
	FILE *stream);

void
sfree(void *p,size_t size);

void
report_memuse(void);

size_t
sfread(void *ptr, size_t size, size_t nitems,  FILE *stream);
