/*

gillStats.h

This is a bunch of statistics code that was given to me by Gill.
Most of it comes from the numerical recipies book and some of it
might be from Nir

*/

#ifndef BEDLONG_H
#define BEDLONG_H

#ifndef BED_H
#include "bed.h"
#endif

#ifndef COMMON_H
#include "common.h"
#endif

struct bedLong
/* Browser extensible data */
{
	struct bedLong *next;   /* Next in singly linked list. */
	char *chrom;	/* Human chromosome or FPC contig */
	long chromStart;	/* Start position in chromosome */
	long chromEnd;	/* End position in chromosome */
	char *name;
	struct slName *goTerms;    /* List of GO Terms */
	char strand;
};

long stringToLong(char *s);

struct bedLong *bedLongLoadN(char *row[], int wordCount);

struct bedLong *filenameToBedLong(char *filename);

struct bedLong *bedToBedLong(struct bed *futon, boolean hasGoTerms);

struct bedLong *cloneBedLong(struct bedLong *futon);

struct bedLong *bedListToBedLong(struct bed *bedList, boolean hasGoTerms);

struct bedLong *cloneBedLongList(struct bedLong *head);

void bedLongFree(struct bedLong **pEl);

void bedLongFreeList(struct bedLong **pList);

void bedLongLineOut(struct bedLong *futon);

void bedLongPrettyOut(struct bedLong *futon);

void showBedLongList(struct bedLong *bedLongList);

struct slName *extractUniqGoTermsFromBedLong(struct bedLong *bedLongList);

boolean bedLongHasGoTerm(struct bedLong *bedLong, char *goTerm);

#endif
