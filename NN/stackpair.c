#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#define INITIAL_SIZE 100

static struct {
	int *ptr;	//Stack top
	int *stack;	//Stack base
	int *end;	//Location after the last stack element
} s_stack;



void stackpairinit_(int *size)
{
	if (s_stack.stack)
		return; //stack already initialized

	s_stack.stack = malloc(2 * INITIAL_SIZE * sizeof(int));

	if (s_stack.stack == 0) {
		perror("Error allocating stack memory. Aborting");
		exit(EXIT_FAILURE);
	}
	
	s_stack.end = s_stack.stack + INITIAL_SIZE;
	s_stack.ptr = s_stack.stack;
}

static void growstack_()
{
	int *tmp = s_stack.stack;
	int oldsize = (s_stack.end - s_stack.stack);
	int newsize = (oldsize + oldsize); // x 2 - totally arbitrary. 

	printf("Growing stack to %d... Consider using a bigger stack right from the start.\n", newsize);
	stackinit_(&newsize);
	memcpy(s_stack.stack, tmp, oldsize * sizeof(int));
	s_stack.ptr += oldsize;
}

void pushpair_(int *obja, int *objb)
{
	assert(s_stack.stack);	

	if (s_stack.ptr == s_stack.end)
		growstack_();		//this is expensive, avoid

	*s_stack.ptr++ = *obja;
	*s_stack.ptr++ = *objb;
}

void poppair_(int *obja, int *objb)
{
	assert(s_stack.stack);
	assert(s_stack.ptr > s_stack.stack);

	*objb = *--s_stack.ptr;
	*obja = *--s_stack.ptr;
}

void stackpairempty_(int *is_empty)
{
	*is_empty = (s_stack.ptr == s_stack.stack) ? 1 : 0;
}

void stackpairflush_()
{
	if (s_stack.ptr == s_stack.stack)
		return;
	printf("flushing stack...\n");
	s_stack.ptr = s_stack.stack;
}

void stackpairdestroy_()
{
	if (s_stack.stack) {
		free(s_stack.stack);
		s_stack.stack = 0;
		s_stack.end = 0;
		s_stack.ptr = 0;
	}
}
