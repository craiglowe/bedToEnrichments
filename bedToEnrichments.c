/*

bedToEnrichments.c
Written by Craig Lowe

*/

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "memalloc.h"
#include "bed.h"
#include "bedLong.h"
#include "dystring.h"
#include "gsl/gsl_cdf.h"


/*---------------------------------------------------------------------------*/

static struct optionSpec optionSpecs[] =
/* command line option specifications */
{
	{"geneAssignments", OPTION_BOOLEAN},
	{"binom", OPTION_BOOLEAN},
	{"hypergeo", OPTION_BOOLEAN},
	{"bonferroni", OPTION_BOOLEAN},
	{"maxExpansion", OPTION_INT},
	{"noExpansionOverlap", OPTION_BOOLEAN},
	{"maxPvalue", OPTION_DOUBLE},
	{"guessTxStart", OPTION_BOOLEAN},
	{"goTermToEnglish", OPTION_STRING},
	{"showNames", OPTION_BOOLEAN},
	{"showParams", OPTION_BOOLEAN},
	{"largeSet", OPTION_STRING},
	{"countUnassigned", OPTION_BOOLEAN},
	{NULL, 0}
};


boolean optGeneAssignments = FALSE;
boolean optBinom = FALSE;
boolean optHypergeo = FALSE;
boolean optBonferroni = FALSE;
int optMaxExpansion = 1000000;
boolean optNoExpansionOverlap = FALSE;
double optMaxPvalue = 0.05;
boolean optGuessTxStart = FALSE;
char *optGoTermToEnglish = NULL;
boolean optShowNames = FALSE;
boolean optShowParams = FALSE;
char *optLargeSet = NULL;
boolean optCountUnassigned = FALSE;


/*---------------------------------------------------------------------------*/

void usage()
/* Explain usage and exit. */
{
errAbort(
	"bedToEnrichments - do enrichment tests when given a .bed file.\n"
	"usage:\n"
	"   bedToEnrichments elements.bed genes.bedLong noGaps.bed\n"
	"options:\n"
	"   -binom                FALSE    use the binomial method\n"
	"   -hypergeo             FALSE    use the hypergeometric method\n"
	"   -bonferroni           FALSE    correct pvalues for multiple tests\n"
	"   -maxExpansion=int     1000000  element will not be assigned to a gene if it is further away than this\n"
	"   -noExpansionOverlap   FALSE    expansion can only happen into bases that have not been assigned to another gene\n"
	"   -maxPvalue=double     0.05     do not print pvalues that are greater than this cutoff\n"
	"   -guessTxStart         FALSE    convert the interval into a point based on the strand information\n"
	"   -goTermToEnglish=str  NULL     file mapping goTerms to english definitions\n"
	"   -showNames            FALSE    show the names of genes hit in the output\n"
	"   -showParams           FALSE    show the parameters used to calculate the p-value\n"
	"   -largeSet=str.bed     NULL     a larger bed file that contains the bases from elements.bed.  This is used like a null model\n"
	"   -geneAssignments      FALSE    just show the elements and the genes assigned to it\n"
	"   -countUnassigned      FALSE    count the elements outside of maxExpansion when doing stats\n"
	"notes:\n"
	"   genes.bedLong is the same format as a 6 column bed, but the score field is replaced with a\n"
	"     comma separated list of GO terms\n"
	"references:\n"
	"  This code has been used and described in:\n"
	"    Lowe CB, Kellis M, Siepel A, Raney BJ, Clamp M, Salama SR, Kingsley DM, Lindblad-Toh K, Haussler D.\n"
	"    Three periods of Regulatory Innovation During Vertebrate Evolution\n"
	"    Science. 2011 Aug 19;333(6045):1019-24.\n"
	"    PMID: 21852499\n"
	"     -and-\n"
	"    Lowe CB, Bejerano G, Haussler D.\n"
	"    Thousands of human mobile element fragments undergo strong purifying selection near developmental genes.\n"
	"    Proc Natl Acad Sci U S A. 2007 May 8;104(19):8005-10.\n"
	"    PMID: 17463089\n"
	);
}


/*---------------------------------------------------------------------------*/


struct slNameDouble
{
	struct slNameDouble *next;
	char* name;
	double number;
};


struct slNameDouble *createSlNameDouble(char *name, double number)
{
	struct slNameDouble *ret = NULL;
	AllocVar(ret);
	ret->name = cloneString(name);
	ret->number = number;
	return(ret);
}


int slNameDoubleCmp(const void *va, const void *vb)
{
	const struct slNameDouble *a = *((struct slNameDouble **)va);
	const struct slNameDouble *b = *((struct slNameDouble **)vb);
	if(a->number > b->number){return(1);}
	else if(a->number == b->number){return(0);}
	else{return(-1);}
}


void showNameDoubleList(struct slNameDouble *head)
{
	struct slNameDouble *curr;
	for(curr=head; curr != NULL; curr=curr->next)
	{
		fprintf(stdout,"%s\t%g\n",curr->name,curr->number);
	}
}


struct hash *fileLoadHash(char *fileName) 
{
	struct hash *hash = newHash(9);
	struct lineFile *lf = lineFileOpen(fileName, TRUE);
	char *row[2];

	while (lineFileRowTab(lf, row))
	{
		hashAdd(hash, row[0], cloneString(row[1]));
	}

	lineFileClose(&lf);
	return hash;
}


char *slNameToCommaString(struct slName *head)
{
	int spaceNeeded = 0;
	struct slName *curr = NULL;
	char *answer = NULL;

	for(curr = head; curr != NULL; curr=curr->next)
	{
		spaceNeeded += strlen(curr->name);
		spaceNeeded += 1;
	}

	AllocArray(answer,spaceNeeded);

	for(curr = head; curr != NULL; curr=curr->next)
	{
		strcat(answer,curr->name);
		if(curr->next != NULL)
			strcat(answer,",");
	}

	return(answer);
}


char *hyperParamsToTabString(long whiteBallsPicked, long totalPicks, long whiteBalls, long totalBalls)
{
	struct dyString *string = newDyString(256);
	char *toRet = NULL;
	double expected = ((double)whiteBalls) / ((double)totalBalls) * ((double)totalPicks);
	dyStringPrintf(string, "%g\t%ld\t%ld\t%ld\t%ld", expected, whiteBallsPicked, totalPicks, whiteBalls, totalBalls);
	toRet = cloneString(string->string);
	dyStringFree(&string);
	return toRet;
}


char *binomParamsToTabString(double prob, long whiteBallsPicked, long totalPicks)
{
	struct dyString *string = newDyString(256);
	char *toRet = NULL;
	double expected = prob * ((double)totalPicks);
	dyStringPrintf(string, "%g\t%ld\t%ld\t%g\tfoo", expected, whiteBallsPicked, totalPicks, prob);
	toRet = cloneString(string->string);
	dyStringFree(&string);
	return toRet;
}


char *hashKeyToString(struct hash *hash, char *key)
{
	struct slName *curr = NULL, *head = NULL;
	struct hashEl *el = NULL;
	char *answer = NULL;

	el = hashLookup(hash, key);
	if(el == NULL)
	{
		//errAbort("%s not found in hash",key);
		return(NULL);
	}
	head = newSlName((char *)el->val);
	el = hashLookupNext(el);

	while(el != NULL)
	{
		curr = newSlName((char *)el->val);
		slAddHead(&head,curr);
		el = hashLookupNext(el);
	}

	slReverse(head);
	answer = slNameToCommaString(head);
	slNameFreeList(&head);
	return(answer);
}


void displayResults(struct slNameDouble *head, struct hash *hitsHash, struct hash *paramsHash)
{
	slSort(&head,slNameDoubleCmp);
	struct hash *goToEnglishHash = NULL;

	if(optGoTermToEnglish != NULL)
		goToEnglishHash = fileLoadHash(optGoTermToEnglish);

	struct slNameDouble *curr;
	for(curr=head; curr != NULL; curr=curr->next)
	{
		if(curr->number <= optMaxPvalue)
		{
		struct dyString *string = newDyString(256);
		dyStringPrintf(string,"%s\t%g",curr->name,curr->number);

		if(paramsHash != NULL){dyStringPrintf(string,"\t%s",(char *)hashMustFindVal(paramsHash,curr->name));}
		if(goToEnglishHash != NULL){dyStringPrintf(string,"\t%s",(char *)hashMustFindVal(goToEnglishHash,curr->name));}
		if(hitsHash != NULL)
		{
			dyStringPrintf(string,"\t%s",hashKeyToString(hitsHash,curr->name));
		}

		fprintf(stdout,"%s\n",string->string);
		dyStringFree(&string);
		}
	}
}


void bonferroniCorrection(struct slNameDouble *head, int numberOfTests)
{
	struct slNameDouble *curr;
	double dubNOT = (double)numberOfTests;
	for(curr=head; curr != NULL; curr=curr->next)
	{
		curr->number *= dubNOT;
		if(curr->number > 1){curr->number = 1;}
	}
}


void bedLongGuessTxStart(struct bedLong *bedLongList)
{
	struct bedLong *futon = NULL;
	for(futon=bedLongList; futon != NULL; futon=futon->next)
	{
		if(futon->strand == '+')
			futon->chromEnd = futon->chromStart + 1;
		else if (futon->strand == '-')
			futon->chromStart = futon->chromEnd - 1;
		else
			errAbort("tried to guess the txStart when there is not strand %s %ld %ld", futon->chrom, futon->chromStart, futon->chromEnd);
	}
}


void expandBedLongListByDistance(struct bedLong *bedLongList, long distance)
{
	struct bedLong *futon = NULL;

	for(futon=bedLongList; futon != NULL; futon=futon->next)
	{
		futon->chromStart = max(0,futon->chromStart - distance);
		futon->chromEnd += distance;
	}
}

void expandBedLongListToNeighbor(struct bedLong *bedLongList, long distance)
{
	struct bedLong *prev = NULL, *curr = NULL;
	long middle = 0;

	for(curr=bedLongList; curr != NULL; curr=curr->next)
	{
		if((prev != NULL) && (strcmp(prev->chrom,curr->chrom) != 0))
		{
			prev->chromEnd += distance;
			prev = NULL;
		}

		if(prev==NULL)
		{
			curr->chromStart = max(0,curr->chromStart - distance);
			prev = curr;
			if(curr->next == NULL){curr->chromEnd += distance;}
		}
		else if(curr->chromStart - prev->chromEnd >= 2 * distance)
		{
			prev->chromEnd += distance;
			curr->chromStart = max(0,curr->chromStart - distance);
			prev = curr;
			if(curr->next == NULL){curr->chromEnd += distance;}
		}
		else if(curr->chromStart - prev->chromEnd >= 0)
		{
			middle = (curr->chromStart + prev->chromEnd)/2;
			prev->chromEnd = middle;
			curr->chromStart = middle;
			prev = curr;
			if(curr->next == NULL){curr->chromEnd += distance;}
		}
		else if(curr->chromEnd - prev->chromEnd >= 0)
		{
			prev = curr;
			if(curr->next == NULL){curr->chromEnd += distance;}
		}
		else if(curr->chromEnd < prev->chromEnd)
		{
			if(curr->next == NULL){curr->chromEnd += distance;}
		}
		else
		{
			errAbort("should not exhaust this if statement");
		}
	}
}


void showSlNameList(struct slName *slNameList)
{
	struct slName *craig;

	for(craig=slNameList; craig != NULL; craig=craig->next)
	{
		fprintf(stdout,"%s\n",craig->name);
	}
}


int bedLongCmp(const void *va, const void *vb)
{
	const struct bedLong *a = *((struct bedLong **)va);
	const struct bedLong *b = *((struct bedLong **)vb);
	int dif;
	dif = strcmp(a->chrom, b->chrom);
	if(dif != 0){return(dif);}
	else if(a->chromStart > b->chromStart){return(1);}
	else if(a->chromStart == b->chromStart){return(0);}
	else{return(-1);}
}


int bedLongCmpStart(struct bedLong *futon, struct bedLong *bunk)
{
	int diff = 0;
	diff = strcmp(futon->chrom, bunk->chrom);
	if(diff == 0)
	{
		if(futon->chromStart < bunk->chromStart){return(-1);}
		else if(futon->chromStart > bunk->chromStart){return(1);}
		else{return(0);}
	}
	else{return(diff);}
}


int bedLongCmpEnd(struct bedLong *futon, struct bedLong *bunk)
{
	int diff = 0;
	diff = strcmp(futon->chrom, bunk->chrom);
	if(diff == 0)
	{
		if(futon->chromEnd < bunk->chromEnd){return(-1);}
		else if(futon->chromEnd > bunk->chromEnd){return(1);}
		else{return(0);}
	}
	else{return(diff);}
}


boolean bedLongOverlap(struct bedLong *futon, struct bedLong *bunk)
{
	assert(futon != NULL);
	if(strcmp(futon->chrom,bunk->chrom) == 0)
	{
		if(min(futon->chromEnd,bunk->chromEnd) - max(futon->chromStart,bunk->chromStart) > 0)
		{
		return(TRUE);
		}
	}
	return(FALSE);
}


long bedLongIntersectGoBases(struct bedLong *geneList, char *goTerm, struct bedLong *allowedRegionsList)
{
	/* returns the number of bases in the intersection of the two bed files */
	/* where only the records of geneList that match "goTerm" will be used */
	/* both the bed lists should be sorted with bedLongCmpStart */
	struct bedLong *gene = NULL, *sequenced = NULL;
	char *prevChr = NULL;
	long sum = 0, prevEnd = 0, overlapStart = 0, overlapEnd = 0;

	gene = geneList;
	sequenced = allowedRegionsList;
	prevChr = cloneString(gene->chrom);

	while(gene != NULL && sequenced != NULL)
	{
		if(bedLongHasGoTerm(gene, goTerm))
		{
			if(strcmp(prevChr,gene->chrom) != 0)
			{
				prevChr = cloneString(gene->chrom);
				prevEnd = 0;
			}

			if(bedLongOverlap(gene,sequenced))
			{
				overlapStart = max(gene->chromStart,sequenced->chromStart);
				overlapEnd = min(gene->chromEnd,sequenced->chromEnd);
				if(overlapStart >= prevEnd)
				{
					sum += overlapEnd - overlapStart;
				}
				else if(overlapEnd > prevEnd)
				{
					sum += overlapEnd - prevEnd;
				}
				prevEnd = max(prevEnd,overlapEnd);
			}
			if(bedLongCmpEnd(gene,sequenced) <= 0){gene = gene->next;}
			else{sequenced = sequenced->next;}
		}
		else
		{
			gene = gene->next;
		}
	}
	return(sum);
}

int bedLongIntersectThreeGoCount(struct bedLong *listOne, char *goTermOne, struct bedLong *listTwo, char *goTermTwo, struct bedLong *listThree, char *goTermThree)
{
	/* returns the number of elements in list one that overlap both something in list */
	/* two and in list three */
	struct bedLong *bedLongOne = NULL, *bedLongTwo = NULL, *bedLongThree = NULL;
	int count = 0;

	bedLongOne = listOne;
	bedLongTwo = listTwo;
	bedLongThree = listThree;

	while(bedLongOne != NULL && bedLongTwo != NULL && bedLongThree != NULL)
	{
		if((goTermOne == NULL) || (bedLongHasGoTerm(bedLongOne, goTermOne)))
		{
			if((goTermTwo == NULL) || (bedLongHasGoTerm(bedLongTwo, goTermTwo)))
			{
				if((goTermThree == NULL) || (bedLongHasGoTerm(bedLongThree, goTermThree)))
				{
					if(bedLongOverlap(bedLongOne,bedLongTwo) && bedLongOverlap(bedLongOne,bedLongThree))
					{
						count++;
						bedLongOne = bedLongOne->next;
					}
					else if((bedLongCmpEnd(bedLongOne,bedLongTwo) < 0) && (bedLongCmpEnd(bedLongOne,bedLongThree) < 0)){bedLongOne = bedLongOne->next;}
					else if(bedLongCmpEnd(bedLongTwo,bedLongThree) < 0){bedLongTwo = bedLongTwo->next;}
					else{bedLongThree = bedLongThree->next;}
				}
				else
				{
					bedLongThree = bedLongThree->next;
				}
			}
			else
			{
				bedLongTwo = bedLongTwo->next;
			}
		}
		else
		{
		bedLongOne = bedLongOne->next;
		}
	}
	return(count);
}

int bedLongIntersectGoCount(struct bedLong *listOne, char *goTermOne, struct bedLong *listTwo, char *goTermTwo, struct hash *retHitsOneHash, struct hash *retHitsTwoHash)
{
	/* returns the number of elements that overlap a gene with the goTerm */
	/* both the bed lists should be sorted with bedLongCmp */
	struct bedLong *bedLongOne = NULL, *bedLongTwo = NULL;
	int count = 0;

	bedLongOne = listOne;
	bedLongTwo = listTwo;

	if((retHitsOneHash != NULL && goTermOne == NULL) || (retHitsTwoHash != NULL && goTermTwo == NULL))
		errAbort("request for names hit, but no go term to use as key");

	while(bedLongOne != NULL && bedLongTwo != NULL)
	{
		if((goTermOne == NULL) || (bedLongHasGoTerm(bedLongOne, goTermOne)))
		{
			if((goTermTwo == NULL) || (bedLongHasGoTerm(bedLongTwo, goTermTwo)))
			{
				if(bedLongOverlap(bedLongOne,bedLongTwo))
				{
					if(retHitsOneHash != NULL)
					{
						if(bedLongOne->name == NULL){errAbort("Error: told to list names, but hit has not name");}
						hashAdd(retHitsOneHash, goTermOne, cloneString(bedLongOne->name));
					}
					if(retHitsTwoHash != NULL)
					{
						if(bedLongTwo->name == NULL){errAbort("Error: told to list names, but hit has not name");}
						hashAdd(retHitsTwoHash, goTermTwo, cloneString(bedLongTwo->name));
					}
					count++;
					bedLongOne = bedLongOne->next;
				}
				else if(bedLongCmpEnd(bedLongOne,bedLongTwo) < 0){bedLongOne = bedLongOne->next;}
				else{bedLongTwo = bedLongTwo->next;}
			}
			else
			{
				bedLongTwo = bedLongTwo->next;
			}
		}
		else
		{
			bedLongOne = bedLongOne->next;
		}
	}
	return(count);
}


int bedLongIntersectCount(struct bedLong *listOne, struct bedLong *listTwo)
{
	/* returns the number of elements from list one that have any overlap with list two */
	/* both the bed lists should be sorted with bedLongCmp */
	struct bedLong *futon = NULL, *bunk = NULL;
	int count = 0;

	futon = listOne;
	bunk = listTwo;

	while(futon != NULL && bunk != NULL)
	{
		if(bedLongOverlap(futon,bunk))
		{
			count++;
			futon = futon->next;
		}
		else if(bedLongCmpEnd(futon,bunk) < 0){futon = futon->next;}
		else{bunk = bunk->next;}
	}
	return(count);
}


long bedLongIntersectBases(struct bedLong *bedLongListA, struct bedLong *bedLongListB)
/* returns the number of bases in the intersection of the two bed files */
/* both the bed lists should be sorted with bedCmp */
{
	struct bedLong *futon = NULL, *bunk = NULL;
	char *prevChr = NULL;
	long sum = 0, prevEnd = 0, overlapStart = 0, overlapEnd = 0;

	futon = bedLongListA;
	bunk = bedLongListB;
	prevChr = cloneString(futon->chrom);

	while(futon != NULL && bunk != NULL)
	{
		if(strcmp(prevChr,futon->chrom) != 0)
		{
			prevChr = cloneString(futon->chrom);
			prevEnd = 0;
		}

		if(bedLongOverlap(futon,bunk))
		{
			overlapStart = max(futon->chromStart,bunk->chromStart);
			overlapEnd = min(futon->chromEnd,bunk->chromEnd);
			if(overlapStart >= prevEnd)
			{
				sum += overlapEnd - overlapStart;
			}
			else if(overlapEnd > prevEnd)
			{
				sum += overlapEnd - prevEnd;
			}
			prevEnd = max(prevEnd,overlapEnd);
		}
		if(bedLongCmpEnd(futon,bunk) <= 0){futon = futon->next;}
		else{bunk = bunk->next;}
	}
	return(sum);
}

/* all poker chips in the bag */
long bedLongBases(struct bedLong *bedLongList)
/* sum the length of all bed files */
/* overlap should only be counted once */
/* bed list should be sorted with bedCmp */
{
	struct bedLong *futon = NULL;
	long sum = 0, prevEnd = 0;
	char *prevChr = NULL;
	prevChr = cloneString(bedLongList->chrom);

	for(futon=bedLongList; futon != NULL; futon=futon->next)
	{
		if(strcmp(prevChr,futon->chrom) != 0)
		{
		prevChr = cloneString(futon->chrom);
		prevEnd = 0;
		}

		if(futon->chromStart > prevEnd)
		{
			sum += (futon->chromEnd - futon->chromStart);
		}
		else if (futon->chromEnd > prevEnd)
		{
			sum += (futon->chromEnd - prevEnd);
		}
		prevEnd = max(prevEnd,futon->chromEnd);
	}
	return(sum);
}


int countGoTermAppearanceInBedLong(struct bedLong *head, char *goTerm)
{
	struct bedLong *gene = NULL;
	int count = 0;
	for(gene=head; gene != NULL; gene=gene->next)
	{
		if(bedLongHasGoTerm(gene, goTerm)){count++;}
	}
	return(count);
}


struct bedLong *findNameInBedLongList(struct bedLong *head, char *name)
{
	struct bedLong *curr = NULL;
	for(curr=head; curr!=NULL; curr=curr->next)
	{
		if(sameString(curr->name, name))
		{
			return(curr);
		}
	}
	return(NULL);
}


long int absDiff(long int a, long int b)
{
	if(a >= b){return(a-b);}
	else{return(b-a);}
}


long int distanceBetweenBeds(struct bedLong *a, struct bedLong *b)
{
	if(differentString(a->chrom, b->chrom)){errAbort("Error: can not calculate distance between beds on separate chroms");}
	if(bedLongOverlap(a,b)){return(0);}
	else
	{
		return(min(absDiff(a->chromStart, b->chromEnd-1), absDiff(a->chromEnd-1, b->chromStart)));
	}
}


struct slNameDouble *hypergeometricNullModelStyle(struct bedLong *elementsList, struct bedLong *largeSetList, struct bedLong *genesList, struct slName *goTerms, struct bedLong *okRegionsList, struct hash *paramsHash)
{
	int totalBalls = 0, whiteBalls = 0, totalPicks = 0, whiteBallsPicked = 0;
	struct slName *term = NULL;
	double pValue = 0;
	struct slNameDouble *termAndPvalue = NULL;

	verbose(2,"  Calculating numbers that will not change in loop\n");
	totalBalls = slCount(largeSetList);
	totalPicks = bedLongIntersectCount(largeSetList,elementsList);

	verbose(2,"  Entering Loop\n");
	for(term=goTerms; term!=NULL; term=term->next)
	{
		whiteBalls = bedLongIntersectGoCount(largeSetList, NULL, genesList, term->name, NULL, NULL);
		whiteBallsPicked = bedLongIntersectThreeGoCount(largeSetList, NULL, genesList, term->name, elementsList, NULL);
		if(paramsHash != NULL){hashAdd(paramsHash,term->name,hyperParamsToTabString(whiteBallsPicked,totalPicks,whiteBalls,totalBalls));}
		//pValue = hyperGeoPValue(whiteBallsPicked, totalPicks, whiteBalls, totalBalls);
		if(whiteBallsPicked == 0){pValue = 1;}
		else{pValue = gsl_cdf_hypergeometric_Q((unsigned int)whiteBallsPicked-1, (unsigned int)whiteBalls, (unsigned int)totalBalls-whiteBalls, (unsigned int)totalPicks);}
		struct slNameDouble *temp = createSlNameDouble(term->name,pValue);
		slAddHead(&termAndPvalue,temp);
	}
	verbose(2,"  Done With Loop\n");

	return(termAndPvalue);
}


struct slNameDouble *hypergeometricStyle(struct bedLong *elementsList, struct bedLong *genesList, struct slName *goTerms, struct bedLong *okRegionsList, struct hash *retHitsHash, struct hash *paramsHash)
{
	int totalBalls = 0, whiteBalls = 0, totalPicks = 0, whiteBallsPicked = 0;
	struct slName *term = NULL;
	double pValue = 0;
	struct slNameDouble *termAndPvalue = NULL;

	verbose(2,"  Calculating numbers that will not change in loop\n");
	totalBalls = slCount(genesList);

	verbose(2,"  Entering Loop\n");
	for(term=goTerms; term!=NULL; term=term->next)
	{
		whiteBalls = countGoTermAppearanceInBedLong(genesList,term->name);
		totalPicks = bedLongIntersectCount(genesList,elementsList);
		whiteBallsPicked = bedLongIntersectGoCount(genesList, term->name, elementsList, NULL, retHitsHash, NULL);
		if(paramsHash != NULL){hashAdd(paramsHash,term->name,hyperParamsToTabString(whiteBallsPicked,totalPicks,whiteBalls,totalBalls));}
		//pValue = hyperGeoPValue(whiteBallsPicked, totalPicks, whiteBalls, totalBalls);
		if(whiteBallsPicked == 0){pValue = 1;}
		else{pValue = gsl_cdf_hypergeometric_Q((unsigned int)whiteBallsPicked-1, (unsigned int)whiteBalls, (unsigned int)totalBalls-whiteBalls, (unsigned int)totalPicks);}
		struct slNameDouble *temp = createSlNameDouble(term->name,pValue);
		slAddHead(&termAndPvalue,temp);
	}
	verbose(2,"  Done With Loop\n");

	return(termAndPvalue);
}


struct slNameDouble *binomialStyle(struct bedLong *elementsList, struct bedLong *genesList, struct slName *goTerms, struct bedLong *okRegionsList, struct hash *retHitsHash, struct hash *paramsHash)
{
	long totalBalls = 0, whiteBalls = 0, totalPicks = 0, whiteBallsPicked = 0;
	struct slName *term = NULL;
	double prob = 0, pValue = 0;
	struct slNameDouble *termAndPvalue = NULL;

	verbose(2,"  Calculating numbers that will not change in loop\n");
	totalBalls = bedLongBases(okRegionsList);
	//totalPicks = slCount(elementsList);
	if(optCountUnassigned){totalPicks = slCount(elementsList);}
	else{totalPicks = bedLongIntersectCount(elementsList,genesList);}

	verbose(2,"  Entering Loop\n");
	for(term=goTerms; term!=NULL; term=term->next)
	{
		whiteBalls = bedLongIntersectGoBases(genesList, term->name, okRegionsList);
		whiteBallsPicked = bedLongIntersectGoCount(elementsList, NULL, genesList, term->name, NULL, retHitsHash);
		prob = ((double)whiteBalls)/((double)totalBalls);
		if(paramsHash != NULL){hashAdd(paramsHash,term->name,binomParamsToTabString(prob,whiteBallsPicked,totalPicks));}
		//pValue = binomPValue(whiteBallsPicked,totalPicks,prob);
		if(whiteBallsPicked == 0){pValue = 1;}
		else{pValue = gsl_cdf_binomial_Q((unsigned int)whiteBallsPicked-1, prob, (unsigned int)totalPicks);}
		struct slNameDouble *temp = createSlNameDouble(term->name,pValue);
		slAddHead(&termAndPvalue,temp);
	}
	verbose(2,"  Done With Loop\n");

	return(termAndPvalue);
}

void assignmentStyle(struct bedLong *elementsList, struct bedLong *genesList, struct bedLong *okRegionsList, struct bedLong *unexpandedGeneList)
{
	struct bedLong *bedLongOne = NULL, *bedLongTwo = NULL;

	bedLongOne = elementsList;
	bedLongTwo = genesList;

	while(bedLongOne != NULL && bedLongTwo != NULL)
	{
		if(bedLongOverlap(bedLongOne,bedLongTwo))
		{
			fprintf(stdout,"%s\t%ld\t%ld\t%s\t%s\t%ld\n",bedLongOne->chrom, bedLongOne->chromStart, bedLongOne->chromEnd, bedLongOne->name, bedLongTwo->name, distanceBetweenBeds(bedLongOne, findNameInBedLongList(unexpandedGeneList, bedLongTwo->name)));
			bedLongOne = bedLongOne->next;
		}
		else if(bedLongCmpEnd(bedLongOne,bedLongTwo) < 0)
		{
			fprintf(stdout,"%s\t%ld\t%ld\t%s\tNONE\tNONE\n",bedLongOne->chrom, bedLongOne->chromStart, bedLongOne->chromEnd, bedLongOne->name);
			bedLongOne = bedLongOne->next;
		}
		else{bedLongTwo = bedLongTwo->next;}
	}
	while(bedLongOne != NULL)
	{
		fprintf(stdout,"%s\t%ld\t%ld\t%s\tNONE\tNONE\n",bedLongOne->chrom, bedLongOne->chromStart, bedLongOne->chromEnd, bedLongOne->name);
		bedLongOne = bedLongOne->next;
	}
}

/*---------------------------------------------------------------------------*/

void bedToGoStats(char *elementsInFile, char *genesInFile, char *noGapInFile)
{
	struct bedLong *nonexpandedList = NULL, *elementsBedLongList = NULL, *genesBedLongList = NULL, *okRegionsBedLongList = NULL, *largeSet = NULL;
	struct slName *goTerms = NULL;
	struct slNameDouble *results = NULL;
	struct hash *hitsHash = NULL, *paramsHash = NULL;

	elementsBedLongList = filenameToBedLong(elementsInFile);
	genesBedLongList = filenameToBedLong(genesInFile);
	okRegionsBedLongList = filenameToBedLong(noGapInFile);
	if(optLargeSet != NULL){largeSet = filenameToBedLong(optLargeSet);}

	if(optGuessTxStart)
		bedLongGuessTxStart(genesBedLongList);

	slSort(&elementsBedLongList, bedLongCmp);
	slSort(&genesBedLongList, bedLongCmp);
	slSort(&okRegionsBedLongList, bedLongCmp);
	if(optLargeSet != NULL){slSort(&largeSet, bedLongCmp);}

	goTerms = extractUniqGoTermsFromBedLong(genesBedLongList);

	//expand gene list
	nonexpandedList = cloneBedLongList(genesBedLongList);
	if(optMaxExpansion != 0)
	{
		verbose(2,"Expanding list\n");
		long maxExp = (long)optMaxExpansion;
		if(optNoExpansionOverlap)
			expandBedLongListToNeighbor(genesBedLongList,maxExp);
		else
			expandBedLongListByDistance(genesBedLongList,maxExp);
	}

	//showBedLongList(genesBedLongList);

	if(optShowNames)
		hitsHash = newHash(9);

	if(optShowParams)
		paramsHash = newHash(9);

	//do math
	verbose(2,"Calculating Stats...\n");

	if(optGeneAssignments)
		assignmentStyle(elementsBedLongList,genesBedLongList,okRegionsBedLongList,nonexpandedList);
	else if(optBinom)
		results = binomialStyle(elementsBedLongList,genesBedLongList,goTerms,okRegionsBedLongList,hitsHash,paramsHash);
	else if(optHypergeo && optLargeSet)
		results = hypergeometricNullModelStyle(elementsBedLongList,largeSet,genesBedLongList,goTerms,okRegionsBedLongList,paramsHash);
	else if(optHypergeo && !optLargeSet)
		results = hypergeometricStyle(elementsBedLongList,genesBedLongList,goTerms,okRegionsBedLongList,hitsHash,paramsHash);
	else
		errAbort("Error: end of if statement should not be reached");

	if(optBonferroni)
	{
		verbose(2,"Correcting Results For Multiple Tests...\n");
		bonferroniCorrection(results,slCount(goTerms));
	}

	verbose(2,"Displaying Results...\n");
	displayResults(results,hitsHash,paramsHash);

	//bedLongFreeList(&elementsBedLongList);
	//bedLongFreeList(&genesBedLongList);
	//bedLongFreeList(&okRegionsBedLongList);
}

/*---------------------------------------------------------------------------*/


int main(int argc, char *argv[])
/* Process command line. */
{
	optionInit(&argc, argv, optionSpecs);
	if (argc != 4)
		usage();

	optGeneAssignments = optionExists("geneAssignments");
	optBinom = optionExists("binom");
	optHypergeo = optionExists("hypergeo");
	optBonferroni = optionExists("bonferroni");
	optMaxExpansion = optionInt("maxExpansion",optMaxExpansion);
	optNoExpansionOverlap = optionExists("noExpansionOverlap");
	optMaxPvalue = optionDouble("maxPvalue",optMaxPvalue);
	optGuessTxStart = optionExists("guessTxStart");
	optGoTermToEnglish = optionVal("goTermToEnglish", NULL);
	optShowNames = optionExists("showNames");
	optShowParams = optionExists("showParams");
	optLargeSet = optionVal("largeSet", NULL);
	optCountUnassigned = optionExists("countUnassigned");
	if (optBinom && optHypergeo)
		errAbort("You can't use both -binom and -hypergeo");
	if (!optBinom && !optHypergeo && !optGeneAssignments)
		errAbort("You must use either -binom or -hypergeo");
	if (optLargeSet && !optHypergeo)
		errAbort("You must use either -hypergeo with -largeSet");
	if (optLargeSet && optShowNames)
		errAbort("You can not use -showNames with -largeSet");

	bedToGoStats(argv[1],argv[2],argv[3]);
	return 0;
}

