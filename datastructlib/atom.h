#ifndef ATOM_INCLUDED
#define ATOM_INCLUDED
#include <stdint.h>

typedef struct atom{
	struct atom *link;
	int32_t len;
	char *str;
}atom;

extern      int32_t Atom_length(const char* str);
extern const char* Atom_new(const char* str, int32_t len);
extern const char* Atom_string(const char* str);
extern const char* Atom_int32_t(long n);

#endif
