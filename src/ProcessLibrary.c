#include "geneid.h"

extern int CYCLES, NORM;
/* */
void LibraryScan1(packExternalInformation* external, packLib** lib, int Strand, long l1, long l2)
{
	int frame, frameStart, frameEnd, iFrame, pos, k;
	float score;
	float normf;

	score = 0;
	normf = (NORM) ? 100*external->nSequences : external->nSequences;

	if(Strand == FORWARD)
	{
		frameStart = 0;
		frameEnd = FRAMES;
	}
	else
	{
		frameStart = FRAMES;
		frameEnd = 2*FRAMES;
	}
	for(frame = frameStart; frame < frameEnd; frame++)
	{
		memset(external->sr[frame], 0, sizeof(float)*((l2-l1)+1));
	}
	//problem with len[0] -> use l2 instead¿?
	for(pos = 0; pos < l2; pos++)
	{
		if(lib[pos] != NULL)
		{
			for(k = 0; k < lib[pos]->nPairs; k++)
			{
				if(lib[pos]->entry[k] != NULL)
				{
					iFrame = lib[pos]->entry[k]->Frame;
					external->sr[iFrame][pos] += lib[pos]->entry[k]->Score/normf;
				}
			}
		}
	}
}

/* */
void LibraryScan2(float** sr, int Strand, long l1, long l2)
{
	//char mess[MAXSTRING];
	int f, frameStart, frameEnd, pos;

	if(Strand == FORWARD)
	{
		frameStart = 0;
		frameEnd = FRAMES;
	}
	else
	{
		frameStart = FRAMES;
		frameEnd = 2*FRAMES;
	}

	for(f = frameStart; f < frameEnd; f++)
	{
		for(pos = l1+1; pos <= l2; pos++)
		{
			sr[f][pos-(l1+1)] += sr[f][pos-l1];
		}
	}
}

/* */
void LibraryScan(packExternalInformation* external, int Strand, long l1, long l2)
{
	int start,end,pairPos;
	//Position where is acumulating scores, and number of acum
	int ne;
	int frame, frameStart, frameEnd;
	float score;

	score = 0;
	if(Strand == FORWARD)
	{
		frameStart = 0;
		frameEnd = FRAMES;
	}
	else
	{
		frameStart = FRAMES;
		frameEnd = 2*FRAMES;
	}

	/* Fill the regions between l1 & l2 with the scores of the library */
	for(frame = frameStart; frame < frameEnd; frame++)//Process by frame
	{
		memset(external->sr[frame], 0, sizeof(float)*((l2-l1)+1));
		ne = external->procEntry[frame]->ne;
		if(ne > 0)
		{
			/* Skip pairs placed before l1 */
			for(start = 0;
				external->procEntry[frame]->ne > 0 ||
				external->procEntry[frame]->Entry[start] == NULL ||
				external->procEntry[frame]->Entry[start]->ResA < l1 ||
				start < ne;
				start++);

			/* Skip pairs placed after l2 */
			for(end = ne-1;
				external->procEntry[frame]->Entry[end] == NULL ||
				external->procEntry[frame]->Entry[end]->ResA > l2;
				end--);

			/* Pairs of residues between boundaries of l1 and l2*/
			pairPos = external->procEntry[frame]->Entry[start]->ResA;

			for(;start < end; start++)
			{
				//Get scores for this pair
				if(external->procEntry[frame]->Entry[start] != NULL)
				{
					if( pairPos != external->procEntry[frame]->Entry[start]->ResA )
					{
						//external->sr[frame][pairPos - l1] /= external->nSequences;
						pairPos = external->procEntry[frame]->Entry[start]->ResA;
					}
					external->sr[frame][pairPos - l1] += external->procEntry[frame]->Entry[start]->Score/100;
				}
			}
		}
	}
}
/* Declares for the first time the static incidence matrix used in the relaxation process */
int** DeclareIncidenceMatrix(long nSequences, long maxLength)
{
	int i;
	int** mat;

	if((mat = malloc(nSequences * sizeof(int*))) == NULL)
		printError("Declaration error: no enough memory for incidence matrix");
	for(i = 0; i < nSequences; i++)
	{
		if((mat[i] = calloc(maxLength, sizeof(int))) == NULL)
			printError("Declaration error: not enough memory for indicende matrix pointers");
	}
	return mat;
}
/* Clean the positions used in by the pair of sequence-residue */
//packLib**
/* void CleanIncidenceMatrix(int** incidence, int s, int r, packLib** lib) */
void CleanIncidenceMatrix(int** incidence, int s, int r, packExternalInformation* external)
{
	int i,sM,rM;

	if(!incidence)
		printError("Cleaning error: Incidence matrix undeclared");

	incidence[s][r] = 0;
	for(i = 0; i < external->library[s][r]->nPairs; i++)
	{
		sM = external->library[s][r]->entry[i]->SeqB;
		rM = external->library[s][r]->entry[i]->ResB;
		incidence[sM][rM] = 0;
	}
}
/* Scores two residues based in the matches shared between them */
//packLib** & nSeq
/* float ScoreResiduePair(packLib** lib, int nSeq, int s1, int r1, int s2, int r2, int frame) */
float ScoreResiduePair(packExternalInformation* external, int s1, int r1, int s2, int r2, int frame)
{
	int i;
	long sM, rM;
	static int **incidence;
	double score = 0;
	double delta = 0;
	double initial_score = 0;
	float inci;
	float listi;
	//char mess[MAXSTRING];

	if(!incidence)
		incidence = DeclareIncidenceMatrix(external->nSequences,external->maxLength);
	if(external->library[s1][r1] == NULL)
		return 0;
	if(external->library[s2][r2] == NULL)
		return 0;
	incidence[s1][r1] = -999;

	for(i = 0; i < external->library[s1][r1]->nPairs; i++)
	{
		//if(external->library[s1][r1]->entry[i]->Frame == frame)
		if(1)
		{
			sM = external->library[s1][r1]->entry[i]->SeqB;
			rM = external->library[s1][r1]->entry[i]->ResB;
			incidence[sM][rM] = external->library[s1][r1]->entry[i]->Score;
		}
	}
	for(i = 0; i < external->library[s2][r2]->nPairs; i++)
	{
		//if(external->library[s2][r2]->entry[i]->Frame == frame)
		if(1)
		{
			sM = external->library[s2][r2]->entry[i]->SeqB;
			rM = external->library[s2][r2]->entry[i]->ResB;

			if(incidence[sM][rM])
			{
				if(incidence[sM][rM] == -999)
				{
					initial_score = external->library[s2][r2]->entry[i]->Score;
					score += external->library[s2][r2]->entry[i]->Score/NORMF;
				}
				else
				{
					inci = (float)incidence[sM][rM]/NORMF;
					listi = external->library[s2][r2]->entry[i]->Score/NORMF;

					delta = (inci>listi)?inci:listi;
					//delta = MAX(external->library[s2][r2]->entry[i]->Score/NORMF,incidence[sM][rM]/NORMF);
					score += delta;
				}
			}
		}
	}
	//clean structure
	CleanIncidenceMatrix(incidence,s1,r1,external);
	score /= external->nSequences;
	score *= NORMF;

	//sprintf(mess,"diff: %f",initial_score - score);
	//printMess(mess);

	return score;
}
/* Process the raw library and try to improve the scores of the more consistent pairs between sequences */
//int RelaxSequence(packExternalInformation* external, int last, int target)
/* */
int RelaxLibrary(packExternalInformation* external, int last)
{
	int s1,r1,s2,r2,k1,k2,i;
	float s;
	int frame;
	//int nentry;
	//char mess[MAXSTRING];
	int result = 0;
	int* updel = NULL;
	float* w = NULL;
	//int size_list;
	int* tmp_list = NULL;
	float* ftmp_list = NULL;
	int fields = 6;
	int last_entry = 0;

	/* Library sequences */
	for(s1 = 0; s1 < external->nSequences; s1++)
	{
		/* Residues per sequence in the library */
		for(r1 = 0; r1 <= external->lastMatch[s1]; r1++)//* for(frame::FRAME)
		{
			//sprintf(mess,"%d %d",s1,r1);
			//printMess(mess);
			/* Check if the (sequence,residue) has an entry*/
			if(external->library[s1][r1] != NULL)
			{
				/* Check all the pairs added */
				for(k1 = 0; k1 < external->library[s1][r1]->nPairs; k1++)//* external->library[s1][frame]->nPairs
				{
					frame = external->library[s1][r1]->entry[k1]->Frame;
					s2 = external->library[s1][r1]->entry[k1]->SeqB;
					r2 = external->library[s1][r1]->entry[k1]->ResB;
					//if(external->library[s1][r1]->entry[k1]->ResA != -777)
					//{
					if(s2 > s1)
					{

						//get the inverse location
						for(k2 = 0; (k2 < external->library[s2][r2]->nPairs) &&
							(external->library[s2][r2]->entry[k2]->SeqB!=s1 ||
							external->library[s2][r2]->entry[k2]->ResB!=r1 ||
							external->library[s2][r2]->entry[k2]->Frame!=frame); k2++);
						//printMess("d1");
						s = ScoreResiduePair(external,s1,r1,s2,r2,frame);
						//external->library[s1][r1]->entry[k1]->ResA = r1;
						//external->library[s2][r2]->entry[k2]->ResA = r1;

						if((tmp_list = realloc(updel, (last_entry+fields)*sizeof(int))) == NULL)
							printError("Updel memory realloc error");
						if((ftmp_list = realloc(w, ((last_entry/fields)+1)*sizeof(float))) == NULL)
							printError("w memory realloc error");

						updel = tmp_list;
						w = ftmp_list;

						updel[last_entry+S1] = s1;
						updel[last_entry+R1] = r1;
						updel[last_entry+K1] = k1;

						updel[last_entry+S2] = s2;
						updel[last_entry+R2] = r2;
						updel[last_entry+K2] = k2;

						//if(s < 1)
						//{
							//Add to the updel list
							//printMess("do del");
							//sprintf(mess,"Pair deletion: %d %d %d %d",s1, r1, s2, r2);
							//printMess(mess);
							//delete entry entry

							//updel[last_entry+W] = -777;
							//w[last_entry/fields] = -777;
							//DeleteEntryLib(external,s1,r1,k1);
							//DeleteEntryLib(external,s2,r2,k2);

							/*if(external->library[s1][r1]->nPairs == 0)
								k1++;*/
						//}
						//else
						//{
							// It prevents access from seq, res, k...
							/* Need to be individualized the process, here is copying all the matches
							 * but you only have to be left with your sequence
							 */
							//updel[last_entry+W] = s;
							w[last_entry/fields] = s;
							/*
							if(last && s>0)
							{
								//add new entry from the updel list!!
								//but when Up&Del the adreesses will change - realloc
								nentry = external->procEntry[frame]->ne;
								nentry++;
								entryLib** e_tmp;
								//external->procEntry[frame]->Entry[k];
								if((e_tmp = realloc(external->procEntry[frame]->Entry, (nentry)*sizeof(entryLib*))) == NULL )
									printMess("Realloc lost pointer");
								else
									external->procEntry[frame]->Entry = e_tmp;

								external->procEntry[frame]->Entry[nentry-1] = external->library[s1][r1]->entry[k1];
								external->procEntry[frame]->ne++;

								//sprintf(mess,"%d %d %d %d",s1,r1,s2,r2);
								//printMess(mess);
								//e_tmp = NULL;
								//left = DeleteEntryLib(external,s1,r1,k1);
								//if(left == 0)
									//free(external->library[s1][r1]->entry);
							}
							*/
							//external->library[s1][r1]->entry[k1]->Score = s;
							//external->library[s2][r2]->entry[k2]->Score = s;
							//k1++;
						//}
						last_entry += fields;
					//}
					//else
						//k1++;
				}
				}
			}
		}
	}
	for(i = 0; i < last_entry; i+=6)
	{
		UpdateAndDeleteEntryLib(external, &updel[i], w[i/fields]);
	}
	tmp_list = NULL;
	free(updel);
	free(w);

	return result;
}

void ScreenLibraryScores(float** sr, int Strand, long l1, long l2)
{
	//char mess[MAXSTRING];
	int f, frameStart, frameEnd, pos;

	if(Strand == FORWARD)
	{
		frameStart = 0;
		frameEnd = FRAMES;
	}
	else
	{
		frameStart = FRAMES;
		frameEnd = 2*FRAMES;
	}

	for(f = frameStart; f < frameEnd; f++)
	{
		for(pos = l1+1; pos <= l2; pos++)
		{
			sr[f][pos-(l1+1)] += sr[f][pos-l1];
		}
	}
}

/* Call to relaxation process and scaning the library to prepare the scores */
void ProcessLibrary(packExternalInformation* external, packLib** lib, int Strand, long l1, long l2)
{
	//char mess[MAXSTRING];
	//struct timespec requestStart;
	//struct timespec requestEnd;
	//double accum;

	//There is no need for any id of target seq?
	LibraryScan(external, Strand, l1, l2);
	LibraryScan2(external->sr, Strand, l1, l2);

	ScreenLibraryScores(external->sr, Strand, l1, l2);
}

/* Call to relaxation process and scaning the library to prepare the scores */
void ProcessLibrary2(packExternalInformation* external, packLib** lib, int Strand, long l1, long l2)
{
	LibraryScan1(external, lib, Strand, l1, l2);
	LibraryScan2(external->sr, Strand, l1, l2);
}

/* Call to relaxation process and scaning the library to prepare the scores */
//void RelaxLibrary(packExternalInformation* external, char* target_seq)
void PreProcessLibrary(packExternalInformation* external)
{
	int i;
	//char mess[MAXSTRING];
	//char fname[30];
	//struct timespec requestStart;
	//struct timespec requestEnd;
	//double accum;

	//testLibrary("initial_lib.tc_lib", external);

	external->procEntry = calloc(2*FRAMES,sizeof(frameEntry*));

	for(i = 0; i < 2*FRAMES; i++)
	{
		if((external->procEntry[i] = calloc(1,sizeof(frameEntry))) == NULL)
			printError("Not enough memory for process entry list");
	}

	for(i = 0; i < CYCLES; i++)
	{
		//Only relax the sequence to predict
		printMess("Relaxing library");
		RelaxLibrary(external, (i==NCYCLES-1)?1:0);
		//sprintf(fname,"rrelaxed_lib_%d.tc_lib",i);
		//testLibrary(fname, external);
	}
	//printMess("Output the library relaxed");
	//testLibrary(external);
	//testLibrary2("relaxed_lib.tc_lib", external);
	//printMess("Done");
}

