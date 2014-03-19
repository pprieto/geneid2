//Read Library
//Select from library
#include "geneid.h"

long ReadLibrary (char *FileName, packExternalInformation* external)
{
	FILE* file;
	char Locusname[LOCUSLENGTH];
	char line[MAXLINE];
	char mess[MAXSTRING];
	char *copyline;
	int read_pairs = 0;
	int kSeq1,kSeq2;
	int seq1, seq2;
	long res1, res2;
	float score;
	int cons, frame;
	int i, a;
	long nPairs, nSeqs, np, len;
	int *myTranslation;
	nPairs = 0;len = 0;

	//printMess("Start reading library");

	if((file=fopen(FileName,"r"))==NULL)
		printError("The library file can not be opened to read");

	fgets(line,MAXLINE,file);
	if(!strcmp(line,"! TC_LIB_FORMAT_01"))
		printError("Library not using TC_LIB_FORMAT correct");
	else{
		fgets(line,MAXLINE,file);
		sscanf(line,"%ld",&nSeqs);
		myTranslation = calloc(nSeqs, sizeof(int));
		//Alloc for nSeqs in library
		external->library = calloc(nSeqs, sizeof(packLib**));
		if((external->lastMatch = calloc(nSeqs, sizeof(long))) == NULL)
			printError("Not enough memory: array of last maches");
	}
	external->len = calloc((int)nSeqs, sizeof(long));
	external->sr = malloc(2*FRAMES*sizeof(float*));
	external->maxLength = 0;
	external->nSequences = nSeqs;
	for(i = 0; i<nSeqs; i++)
	{
		fgets(line,MAXLINE,file);
		sscanf(line,"%s",Locusname);
		len = strlen(line);
		a = setkeyDict(external->locusNames,Locusname);
		myTranslation[i] = a;
		if(a>MAXSEQLIB)
			printError("Maximum number of sequences in the library exceeded");
		while(line[strlen(line)-1] != '\n')
		{
			fgets(line,MAXLINE,file);
			len += strlen(line);
		}
		//For each seq alloc length;
		external->library[i] = calloc(len, sizeof(packLib*));
		external->len[i] = len;
		if(len > external->maxLength)
			external->maxLength = len;

	}
	for(i = 0; i < 2*FRAMES; i++)
	{
	  external->sr[i] = calloc(external->maxLength, sizeof(float));
	}
	while((fgets(line,MAXLINE,file))!=NULL)
	{
		if(line[0]=='#')
		{
			read_pairs = 1;
			copyline = strtok(line,"#");
			if (sscanf(copyline,"%d %d", &seq1, &seq2) != 2)
				printError("Error reading #sequences in library");
			//sprintf(mess,"#%d %ld",seq1,external->library[kSeq1][res1]->entry[np]->SeqB);
			//printMess(line);
		}
		else if(line[0]=='!')
			printMess("Skipping ending lines");
		else if(read_pairs==1)
		{
			//printMess(line);
			if (sscanf(line,"%ld %ld %f %d %d",&res1, &res2, &score, &cons ,&frame) != 5)
				printError("Error reading first residue of a pair in library");

			kSeq1 = myTranslation[seq1-1];
			kSeq2 = myTranslation[seq2-1];

			if(external->library[kSeq1][res1] == NULL)
			{
				external->library[kSeq1][res1] = calloc(1,sizeof(packLib));
				external->library[kSeq1][res1]->entry = calloc(MAXPAIR,sizeof(entryLib*));
			}

			np = external->library[kSeq1][res1]->nPairs;
			external->library[kSeq1][res1]->entry[np] = RequestNewEntryLib();
			//external->library[kSeq1][res1]->entry[np] = calloc(1,sizeof(entryLib));
			external->library[kSeq1][res1]->entry[np]->ResA = res1;
			external->library[kSeq1][res1]->entry[np]->ResB = res2;
			external->library[kSeq1][res1]->entry[np]->SeqB = kSeq2;
			//Transform the frame from blast to geneid format
			frame = (frame>0)? frame-1 : 3 + (((frame+1)%-3)*-1);

			external->library[kSeq1][res1]->entry[np]->Frame = frame;

			//Normalize the score from the library
			//score = score/external->maxLength;
			external->library[kSeq1][res1]->entry[np]->Score = score;
			//external->library[kSeq1][res1]->entry[np]->Score = score/100;
			external->library[kSeq1][res1]->nPairs++;
			external->library[kSeq1][res1]->avgScore++;
			if(res1 > external->lastMatch[kSeq1])
				external->lastMatch[kSeq1] = res1;

			//sprintf(mess,"#%d %ld",seq1,external->library[kSeq1][res1]->entry[np]->SeqB);
			//printMess(mess);
			sprintf(mess,"%ld %ld %f %d",res1,external->library[kSeq1][res1]->entry[np]->ResB,external->library[kSeq1][res1]->entry[np]->Score, frame);
			//printMess(mess);

			if(external->library[kSeq2][res2] == NULL)
			{
				external->library[kSeq2][res2] = calloc(1,sizeof(packLib));
				external->library[kSeq2][res2]->entry = calloc(MAXPAIR,sizeof(entryLib*));
			}

			np = external->library[kSeq2][res2]->nPairs;
			external->library[kSeq2][res2]->entry[np] = RequestNewEntryLib();
			external->library[kSeq2][res2]->entry[np]->ResA = res2;
			external->library[kSeq2][res2]->entry[np]->ResB = res1;
			external->library[kSeq2][res2]->entry[np]->SeqB = kSeq1;
			external->library[kSeq2][res2]->entry[np]->Frame = frame;

			//Normalize the score from the library
			//external->library[kSeq2][res2]->entry[np]->Score = score;
			//external->library[kSeq2][res2]->entry[np]->Score = score/100;
			external->library[kSeq2][res2]->entry[np]->Score = score;
			external->library[kSeq2][res2]->nPairs++;
			external->library[kSeq2][res2]->avgScore++;
			if(res2 > external->lastMatch[kSeq2])
				external->lastMatch[kSeq2] = res2;

			nPairs++;
		}
		//else printMess("Skipping header lines");
	}
	for(seq1=0; seq1<nSeqs; seq1++)
	{
		kSeq1 = myTranslation[seq1];
		np = 0;
		for(res1=0; res1<external->len[kSeq1]; res1++)
		{
			if(external->library[kSeq1][res1] != NULL && np < external->library[kSeq1][res1]->nPairs)
				external->library[kSeq1][res1]->avgScore/=external->library[kSeq1][res1]->nPairs;
		}
	}
	free(myTranslation);
	return nPairs;
}

void testLibrary(char *FileName, packExternalInformation* external)
{
	FILE* file;
	char mess[MAXSTRING];
	int s1,s2,s2c,np,match;
	long r1, r2, w;
	s2 = 2;

	if((file=fopen(FileName,"w"))==NULL)
		printError("The library file can not be opened to read");

	for(s1 = 0; s1<external->nSequences; s1++)
	{
		sprintf(mess,"SEQUENCE %d\n",s1);
		fprintf(file,mess);

		for(r1 = 0; r1 <= external->lastMatch[s1]; r1++)
		{
			if(external->library[s1][r1] != NULL)
			{
				np = external->library[s1][r1]->nPairs;
			for(match = 0; match<np; match++)
			{
				s2c = (external->library[s1][r1]->entry[match]->SeqB);
					r2 = (long) external->library[s1][r1]->entry[match]->ResB;
					w = (long) external->library[s1][r1]->entry[match]->Score;
					sprintf (mess, "\t%sS1:%5d - R1:%5ld S2:%5d R2:%5ld W:%5ld\n", "", s1, r1, s2c, r2, w);
					fprintf(file,mess);
			}
			}
		}

		}
}

void testLibrary2 (char *FileName, packExternalInformation* external)
{
	FILE* file;
	char mess[MAXSTRING];
	int s1,s2,s2c,np,match;
	long r1, r2, w;
	s2 = 2;

	if((file=fopen(FileName,"w"))==NULL)
		printError("The library file can not be opened to read");

	for(s1 = 0; s1<external->nSequences; s1++)
	{
		for(s2 = s1+1; s2< external->nSequences; s2++)
		{
			sprintf(mess,"#%d %d\n",s1+1,s2+1);
			fprintf(file,mess);
			for(r1 = 0; r1 <= external->lastMatch[s1]; r1++)
			{
				if(external->library[s1][r1] != NULL)
				{
					np = external->library[s1][r1]->nPairs;
				for(match = 0; match<np; match++)
				{
					if(external->library[s1][r1]->entry[match] != NULL)
					{
						s2c = (external->library[s1][r1]->entry[match]->SeqB)+1;
						if(s2+1==s2c)
						{
							r2 = (long) external->library[s1][r1]->entry[match]->ResB;
							w = (long) external->library[s1][r1]->entry[match]->Score;
							sprintf(mess,"%5ld %5ld %5f %5d %5d\n",r1,external->library[s1][r1]->entry[match]->ResB,external->library[s1][r1]->entry[match]->Score,0,external->library[s1][r1]->entry[match]->Frame);
							fprintf(file,mess);
							/*if(external->library[s1][r1]->entry[match]->ResA != -777)
							{
								sprintf(mess, "No -777 %d %ld",s1,r1);
								printMess(mess);
								exit(1);
							}*/
						}
					}
				}
				}
			}

		}
	}
}

packLib** SelectLibrary(packExternalInformation* external, char* Locus, long LengthSequence)
{
	packLib** pLib;
	int a;

	a = getkeyDict(external->locusNames, Locus);
	if (a == NOTFOUND)
		pLib = NULL;
	pLib = external->library[a];

	return pLib;
}

int TranslateLocusToSeqLib(packExternalInformation* external, char* Locus, long LengthSequence)
{
	int seq;

	seq = getkeyDict(external->locusNames, Locus);
	return seq;
}
void UpdateAndDeleteEntryLib(packExternalInformation* external, int* pair_target, float w)
{
	if(w > 1)
	{
		UpdatePairEntryLib(external, pair_target, w);
		//external->library[pair_target[S1]][pair_target[R1]]->entry[pair_target[K1]]->Score = w;
		//external->library[pair_target[S2]][pair_target[R2]]->entry[pair_target[K2]]->Score = w;
	}
	else
	{
		//DeleteEntryLib(external, pair_target[S1], pair_target[R1], pair_target[K1]);
		//DeleteEntryLib(external, pair_target[S2], pair_target[R2], pair_target[K2]);
		DeletePairEntryLib(external, pair_target);
	}
}
void UpdatePairEntryLib(packExternalInformation* external, int* pair_target, float w)
{
	int s1,s2,r1,r2;
	int i;

	s1 = pair_target[S1];
	s2 = pair_target[S2];
	r1 = pair_target[R1];
	r2 = pair_target[R2];

	for(i = 0; i < external->library[s1][r1]->nPairs &&
				(external->library[s1][r1]->entry[i] == NULL ||
				external->library[s1][r1]->entry[i]->SeqB != s2 ||
				external->library[s1][r1]->entry[i]->ResB != r2); i++);

	external->library[s1][r1]->entry[i]->Score = w;

	for(i = 0; i < external->library[s2][r2]->nPairs &&
				(external->library[s2][r2]->entry[i] == NULL ||
				external->library[s2][r2]->entry[i]->SeqB != s1 ||
				external->library[s2][r2]->entry[i]->ResB != r1); i++);

	external->library[s2][r2]->entry[i]->Score = w;
}


void DeletePairEntryLib(packExternalInformation* external, int* pair_target)
{
	int s1,s2,r1,r2;
	int i;

	s1 = pair_target[S1];
	s2 = pair_target[S2];
	r1 = pair_target[R1];
	r2 = pair_target[R2];

	for(i = 0; i < external->library[s1][r1]->nPairs &&
				(external->library[s1][r1]->entry[i] == NULL ||
				external->library[s1][r1]->entry[i]->SeqB != s2 ||
				external->library[s1][r1]->entry[i]->ResB != r2); i++);

	if(DeleteEntryLib(external, s1, r1, i) == 0)
		free(external->library[s1][r1]->entry);

	for(i = 0; i < external->library[s2][r2]->nPairs &&
				(external->library[s2][r2]->entry[i] == NULL ||
				external->library[s2][r2]->entry[i]->SeqB != s1 ||
				external->library[s2][r2]->entry[i]->ResB != r1); i++);

	 if(DeleteEntryLib(external, s2, r2, i) == 0)
		free(external->library[s2][r2]->entry);

}

int DeleteEntryLib(packExternalInformation* external, int s,int r,int k)
{
	int avgm,i,left, initial_np;
	entryLib* e;
	entryLib* del;
	//char mess[MAXSTRING];

	//sprintf(mess,"s:%d r:%d k:%d", s, r, k);
	//printMess(mess);

	//delete the target
	initial_np = external->library[s][r]->nPairs;
	del = external->library[s][r]->entry[k];

	external->library[s][r]->entry[k] = NULL;
	external->library[s][r]->nPairs--;
	left = external->library[s][r]->nPairs;
	/*
	 if(last==0) there is nothing to shift
	 decrement pairs, and realloc only
	*/
	for(i = k ;i < left; i++)
	{
		//shift behind
		e = external->library[s][r]->entry[i+1];
		external->library[s][r]->entry[i] = e;
	}

	if(i > 0 && i != k)
	{
		external->library[s][r]->entry[left] = NULL;
	}
	free(del);
	del = NULL;
	//printMess("d");
	external->library[s][r]->entry = realloc(external->library[s][r]->entry, left*sizeof(entryLib*));

	//update score related stuff
	avgm = (external->library[s][r]->nPairs)?external->library[s][r]->avgScore/external->library[s][r]->nPairs:0;
	external->library[s][r]->avgScore -= avgm;

	return left;
}
void cleanLibraryMemory(packExternalInformation* external)
{
	//free(external->lastMatch);
}
