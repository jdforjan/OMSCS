#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include <math.h>

//From office hours
//head should probably be treated as first sublist head
//calculate sums on each sublist, calculate how far heads are from each other, add head value to rest of sublists

void parallelListRanks(long head, const long* ListSuccessor, long* ListPrefixData, size_t n)
{
	int p = omp_get_num_procs();//get number of processors on this system>
	long s = p*ceil(log(n) / log(2));; //optimal number of sublists
	
	long* sublistHead = (long*)malloc(s * sizeof(long));
	long* sublistPrefixData = (long*)malloc(s * sizeof(long));//array to save all the previous subarray's size
	long* sublistSuccessor = (long*)malloc(s * sizeof(long)); //stores the next head of each subhead, so wo can know the sublist orders
	
	//printf("optimal sublists %ld ", s);
	generateHeadNodes(head, sublistHead, s, n);
	computeListRanking(ListPrefixData, ListSuccessor, sublistHead, sublistPrefixData, sublistSuccessor, s);
	calculatePrefixSums(sublistHead, sublistPrefixData, sublistSuccessor, head, ListPrefixData, s);
	finalizeListRanking(ListPrefixData, ListSuccessor, sublistHead, sublistPrefixData, sublistSuccessor);

	free(sublistHead);
	free(sublistPrefixData);
	free(sublistSuccessor);
}

//DONE
//1. Sequentially, partition the input list into s sublists, by randomly choosing s sublist head nodes.
void generateHeadNodes(long head, long* sublistHead, size_t s, long n)
{
	long randHead;
	long memoryBlock = n / (s - 1);
	//printf("\n head: %ld", head);
	//printf("\n n: %ld", s - 1);
	//printf("\n s: %ld \n", n);
	//printf("\n memoryBlock: %ld \n", memoryBlock);

	//sublistHead[0] = head; 
	for (long i = 0; i < s - 1; i++)//(long i = 1; i < s; i++)
	{
		randHead = (memoryBlock * i) + (rand() % memoryBlock);
		if (randHead == head || randHead >= n) //ensure random head is valid
		{
			i--;
			continue; //restart the i loop
		}
		else
		{
			sublistHead[i] = randHead;
			//printf(" %ld \n", randHead);
		}
	}

	sublistHead[s - 1] = head;  //The last element in sublistHead is the global head
}

//2. Processor pi traverses each sublist, computing the local (with respect to the start of the sublist) ranks of each element and store in R[i].
//2. In parallel, traverse each sublist computing the list ranking of each node within the sublists. AKA .. Form prefix sums in each sublist
void computeListRanking(long* ListPrefixData, const long* ListSuccessor, long* sublistHead, long* sublistPrefixData, long* sublistSuccessor, long s)
{
#pragma omp parallel
	{
	#pragma omp for
		for (long pi = 0; pi < s; pi++)
		{
			ListPrefixData[sublistHead[pi]] = -1; //Change sublist heads to -1
			sublistSuccessor[pi] = -1; //I don't think this is necessary
		}

		#pragma omp for
		for (long pi = 0; pi < s; pi++)
		{
			long successor = ListSuccessor[sublistHead[pi]]; //Get successor for head i
			long rank = 1; //Initialize rank for each processor
			while (successor != -1) //Traverse until the end of the List is reached
			{
				if (ListPrefixData[successor] == -1)//Traverse until next head is reached
				{
					//Do we need this for loop
					for (int j = 0; j < s; j++)
					{
						if (successor == sublistHead[j])//Find out which head node we are at
						{
							sublistSuccessor[pi] = j;
							sublistPrefixData[j] = rank;
							break;
						}
					}
					break;
				}
				ListPrefixData[successor] = rank++;
				successor = ListSuccessor[successor]; //Updates next successor
			}
		}
	}
}

//3. Sequentially, compute the list ranking of the head nodes overall(by doing a prefix sum of the head node ranks). AKA .. calculate the prefix sum of the subsize array
void calculatePrefixSums(long* sublistHead, long* sublistPrefixData, long* sublistSuccessor, long head, long* ListPrefixData, long s)
{
	long sublistIndex = s - 1; //This is the index of the global head
	sublistPrefixData[sublistIndex] = 0;
	ListPrefixData[head] = 0;//Prefix of head seems to always be zero
	//ListPrefixData[sublistHead[sublistIndex]] = 0;//Prefix of head seems to always be zero //Can I substitute this because s -1 is global head location?
	while (sublistSuccessor[sublistIndex] != -1)
	{
		sublistPrefixData[sublistSuccessor[sublistIndex]] = sublistPrefixData[sublistSuccessor[sublistIndex]] + sublistPrefixData[sublistIndex];
		sublistIndex = sublistSuccessor[sublistIndex];
		ListPrefixData[sublistHead[sublistIndex]] = sublistPrefixData[sublistIndex];  //write back the sublist head's rank
	}

	//https://www.cs.cmu.edu/~guyb/papers/Ble93.pdf pg 38
	//void allprefixsums(ListPrefixData, sublistPrefixData, sublistSuccessor)
	//	long sublistIndex = 0
	//	long sum = sublistPrefixData[sublistIndex]
	//	ListPrefixData[sublistIndex] = sum
	//	while (sublistSuccessor[sublistIndex] != -1)
	//	{
	//		sublistIndex = sublistSuccessor[sublistIndex]
	//		sum = sum + sublistPrefixData[sublistIndex]
	//		ListPrefixData[sublistIndex] = sum
	//	}
}

//4. In parallel, traverse the sublists again to convert the sublist ranking to the complete list ranking(by adding the previously determined prefix sum values).
//parallel add back each subhead's rank to each elements in the sublist's rank
void finalizeListRanking(long* ListPrefixData, const long* ListSuccessor, long* sublistHead, long* sublistPrefixData, long* sublistSuccessor)
{
	#pragma omp parallel
	{
		#pragma omp for
		for (long pi = 0; pi < s; pi++)
		{
			long successor = ListSuccessor[sublistHead[pi]]; //Get successor for head i
			while (successor != -1) //Traverse until the end of the sublist is reached
			{
				if (successor == sublistHead[sublistSuccessor[pi]]) //Is this necessary?
				{
					break;
				}
				ListPrefixData[successor] = ListPrefixData[successor] + sublistPrefixData[pi]; //Add the sublist rank to each item in the sublist
				successor = ListSuccessor[successor]; //Updates next successor
			}

			//while ((nextSuccessor != -1) && (nextSuccessor != sublistHead[sublistSuccessor[pi]])) //Traverse until the end of the sublist is reached
			//{
			//	ListPrefixData[nextSuccessor] = ListPrefixData[nextSuccessor] + sublistPrefixData[pi];
			//	nextSuccessor = ListSuccessor[nextSuccessor]; //Updates next successor
			//}
		}
	}
}
