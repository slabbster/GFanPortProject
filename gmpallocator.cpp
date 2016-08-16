#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cstdio> /* Always include cstdio before gmp.h.*/
#include <gmp.h>
#include <assert.h>

#define NBUCKETSGMP 9

/*static inline int bufNum(int size)
{
  if(size<=8)return 0;
  else if(size<=16)return 1;
  else if(size<=32)return 2;
  else if(size<=64)return 3;
  else if(size<=128)return 4;
  else if(size<=256)return 5;
  else if(size<=512)return 6;
  else if(size<=1024)return 7;
  else if(size<=2048)return 8;
  return -1;
}
*/
static inline int bufNum(int size)
{
  if(size<=8)return 0;
  int ret=0;
  while(size>8)
    {
      ret++;
      if(ret>=(NBUCKETSGMP))return -1;
      size>>1;
    }
  return ret;
}

class AllocBucket
{
  void *linkList;
  int size;
public:
  void init(int size_)
  {
    linkList=0;
    size=size_;
  }
  void grow()
  {
    int bufSize=1024*64;
    void *buf=malloc(bufSize);
    assert(buf);
    int n=bufSize/size;
    for(int i=n-1;i>=0;i--)
      {
	void **addr=(void**)((char*)buf+i*size);
	*addr=linkList;
	linkList=addr;
      }
    //    fprintf(stderr,"GROWING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  }
  void *alloc()
  {
    do
      {
	if(linkList)
	  {
	    void *ret=linkList;
	    linkList=*((void**)linkList);
	    return ret;
	  }
	grow();
      }
    while(1);
    return 0;
  }
  void free(void* ptr)
  {
    *((void **)ptr)=linkList;
    linkList=ptr;
  }
};

static AllocBucket theBucketList[NBUCKETSGMP];

void *myAllocateFunction(size_t alloc_size)
{
  //  fprintf(stderr,"Allocating: %i bytes.\n",alloc_size);

  int bn=bufNum(alloc_size);
  if(bn!=-1)
    {
      return theBucketList[bn].alloc();
    }
  return malloc(alloc_size);
}

void *myReallocateFunction (void *ptr, size_t old_size, size_t new_size)
{
  //  fprintf(stderr,"Reallocating: %i --> %i bytes.\n",old_size,new_size);
  if(new_size<=8)return ptr;
  int min=old_size;
  if(new_size<old_size)min=new_size;
  int bn_old=bufNum(old_size);
  int bn_new=bufNum(new_size);
  if(bn_old != -1)
    {
      if(bn_new!=-1)
	{
	  if(bn_new==bn_old)return ptr;
	  else
	    {
	      void *ret=theBucketList[bn_new].alloc();
	      memcpy(ret,ptr,min);
	      theBucketList[bn_old].free(ptr);
	    }
	}
      else
	{
	  void *ret=malloc(new_size);
	  memcpy(ret,ptr,min);
	  theBucketList[bn_old].free(ptr);
	  return ret;
	}
    }
  else
    if(bn_new!=-1)
      {
	void *ret=theBucketList[bn_new].alloc();
	memcpy(ret,ptr,min);
	free(ptr);
	return ret;
      }
  return realloc(ptr,new_size);
}

void myDeallocateFunction (void *ptr, size_t size)
{
  //  fprintf(stderr,"Freeing: %i bytes.\n",size);

  int bn=bufNum(size);
  if(bn!=-1)return theBucketList[bn].free(ptr);
  return free(ptr);
}

//Debug
void *myAllocateFunctionD(size_t alloc_size)
{
  fprintf(stderr,"Allocating: %i bytes.\n",alloc_size);
  return malloc(alloc_size);
}

void *myReallocateFunctionD(void *ptr, size_t old_size, size_t new_size)
{
  fprintf(stderr,"Reallocating: %i --> %i bytes.\n",old_size,new_size);
  return realloc(ptr,new_size);
}

void myDeallocateFunctionD(void *ptr, size_t size)
{
  fprintf(stderr,"Freeing: %i bytes.\n",size);
  return free(ptr);
}


void changeGmpAllocator()
{
  int e=8;
  for(int i=0;i<NBUCKETSGMP;i++)
    {
      theBucketList[i].init(e);
      e+=e;
    }

  mp_set_memory_functions (&myAllocateFunction,&myReallocateFunction,&myDeallocateFunction);
}

class AllocatorDummy
{
public:
  AllocatorDummy()
  {
    changeGmpAllocator();
  }
};

static AllocatorDummy d;//disable allocator for better compatibility
