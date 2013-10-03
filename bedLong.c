/*

bedLong.c

This is a modified version of bed.c from the kent source tree that
was originally started by Jim Kent and has continued to be maintained
and expanded by Jim and others at UCSC.

*/

#include "common.h"
#include "bedLong.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "memalloc.h"
#include "bed.h"


long stringToLong(char *s)
{
	long res = 0;
	char *p, *p0 = s;

	if (*p0 == '-')
		p0++;
	p = p0;
	while ((*p >= '0') && (*p <= '9'))
	{
		res *= 10;
		res += *p - '0';
		p++;
	}
	if ((*p != '\0') || (p == p0))
		errAbort("invalid signed number: \"%s\"", s);
	if (*s == '-')
		return -res;
	else
		return res;
}


struct bedLong *bedLongLoadN(char *row[], int wordCount)
/* Convert a row of strings to a bed. */
{
	struct bedLong *futon;

	AllocVar(futon);
	futon->chrom = cloneString(row[0]);
	futon->chromStart = stringToLong(row[1]);
	futon->chromEnd = stringToLong(row[2]);
	if (wordCount > 3)
		futon->name = cloneString(row[3]);
	if (wordCount > 4)
		futon->goTerms = slNameListFromComma(row[4]);
	if (wordCount > 5)
		futon->strand = row[5][0];
	return futon;
}


struct bedLong *filenameToBedLong(char *filename)
{
	struct bedLong *list = NULL, *el;
	int numFields;
	char *line;
	struct lineFile *lf = lineFileOpen(filename, TRUE);

	lineFileNextReal(lf, &line);
	numFields = countChars(line,'\t') + 1;
	if(numFields < 3 || numFields > 6){errAbort("file %s has %d fields when it needs between 3 and 6",filename,numFields);}
	lineFileClose(&lf);

	lf = lineFileOpen(filename, TRUE);
	char *row[numFields];
	while (lineFileRow(lf, row))
	{
		el = bedLongLoadN(row, numFields);
		slAddHead(&list, el);
	}
	lineFileClose(&lf);
	slReverse(&list);
	return list;
}


struct bedLong *bedToBedLong(struct bed *futon, boolean hasGoTerms)
{
	struct bedLong *ret;
	AllocVar(ret);
	ret->chrom = cloneString(futon->chrom);
	ret->chromStart = futon->chromStart;
	ret->chromEnd = futon->chromEnd;
	if(hasGoTerms)
	{
		ret->name = NULL;
		ret->goTerms = slNameListFromComma(futon->name);
	}
	else
	{
		ret->name = cloneString(futon->name);
		ret->goTerms = NULL;
	}
	return ret;
}


struct bedLong *cloneBedLong(struct bedLong *futon)
{
	struct bedLong *ret;
	AllocVar(ret);
	ret->next=NULL;
	ret->chrom = cloneString(futon->chrom);
	ret->chromStart = futon->chromStart;
	ret->chromEnd = futon->chromEnd;
	ret->name = cloneString(futon->name);
	ret->goTerms = slNameCloneList(futon->goTerms);
	return ret;
}


struct bedLong *bedListToBedLong(struct bed *bedList, boolean hasGoTerms)
{
	struct bedLong *bedLongListOut = NULL;
	struct bed *futon = NULL;

	for (futon=bedList;  futon != NULL;  futon=futon->next)
	{
		struct bedLong *newBedLong = bedToBedLong(futon, hasGoTerms);
		slAddHead(&bedLongListOut, newBedLong);
	}

	slReverse(&bedLongListOut);
	return bedLongListOut;
}


struct bedLong *cloneBedLongList(struct bedLong *head)
{
	struct bedLong *answer = NULL, *curr = NULL;

	for (curr=head; curr != NULL;  curr=curr->next)
	{
		struct bedLong *newBedLong = cloneBedLong(curr);
		slAddHead(&answer, newBedLong);
	}

	slReverse(&answer);
	return answer;
}


void bedLongFree(struct bedLong **pEl)
/* Free a single dynamically allocated bed such as created
 * with bedLoad(). */
{
	struct bedLong *el;

	if ((el = *pEl) == NULL) return;
	freeMem(el->chrom);
	freeMem(el->name);
	slNameFree(el->goTerms);
	freez(pEl);
}


void bedLongFreeList(struct bedLong **pList)
/* Free a list of dynamically allocated bed's */
{
	struct bedLong *el, *next;

	for (el = *pList; el != NULL; el = next)
	{
		next = el->next;
		bedLongFree(&el);
	}
	*pList = NULL;
}


void bedLongLineOut(struct bedLong *futon)
{
	if(futon == NULL){fprintf(stdout,"NULL\n");}
	else if(futon->name == NULL){fprintf(stdout,"%s\t%ld\t%ld\n",futon->chrom,futon->chromStart,futon->chromEnd);}
	else{fprintf(stdout,"%s\t%ld\t%ld\t%s\n",futon->chrom,futon->chromStart,futon->chromEnd,futon->name);}
}


void bedLongPrettyOut(struct bedLong *futon)
{
	if(futon == NULL){fprintf(stdout,"bedLong = NULL\n");}
	else
	{
		if(futon->next == NULL){fprintf(stdout,"  next = NULL\n");}
		if(futon->chrom == NULL){fprintf(stdout,"  chrom = NULL\n");}
		else{fprintf(stdout,"  chrom = %s\n", futon->chrom);}
		fprintf(stdout,"  chromStart = %ld\n", futon->chromStart);
		fprintf(stdout,"  chromEnd = %ld\n", futon->chromEnd);
		if(futon->name == NULL){fprintf(stdout,"  name = NULL\n");}
		else{fprintf(stdout,"  name = %s\n", futon->name);}
	}
}


void showBedLongList(struct bedLong *bedLongList)
{
	struct bedLong *futon = NULL;

	for(futon=bedLongList; futon != NULL; futon=futon->next)
	{
		bedLongLineOut(futon);
	}
}


struct slName *extractUniqGoTermsFromBedLong(struct bedLong *bedLongList)
{
	struct slName *allTerms = NULL, *prevTerms = NULL, *currTerms = NULL;
	struct bedLong *futon = NULL;

	for(futon=bedLongList; futon != NULL; futon=futon->next)
	{
		currTerms = slNameCloneList(futon->goTerms);
		prevTerms = allTerms;
		allTerms = slCat(currTerms,prevTerms);
	}

	slUniqify(&allTerms,slNameCmp,slNameFree);

	return(allTerms);
}


boolean bedLongHasGoTerm(struct bedLong *bedLong, char *goTerm)
{
	return(slNameInList(bedLong->goTerms, goTerm));
}

