/*
 * Tree search generalized from Knuth (6.2.2) Algorithm T just like
 * the AT&T man page says.
 *
 * The node_t structure is for internal use only, lint doesn't grok it.
 *
 * Written by reading the System V Interface Definition, not the code.
 *
 * Totally public domain.
 */
/*LINTLIBRARY*/

#ifndef TSEARCH_H
#define TSEARCH_H

void * tsearch(const void *vkey, void **vrootp,
    int (*compar)(const void *, const void *));

void * tfind(const void *vkey, void **vrootp,
    int (*compar)(const void *, const void *));

void * tdelete(const void *vkey, void **vrootp,
    int (*compar)(const void *, const void *));


#endif