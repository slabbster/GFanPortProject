#include "timer.h"

#include "iostream"
#include <time.h>
#include <assert.h>

using namespace std;

Timer* Timer::timers;

//--------------------------------------------------
// Timer
//--------------------------------------------------

bool Timer::doTimingThisTime()
{
  return false; //Disabling timer
  if(n>=skip)n=0;
  return n==0;
}


Timer::Timer(const char* s, int skip):
  ons(0),
  n(0),
  total(0)
{
  this->s=s;
  this->skip=skip;
  next=timers;
  timers=this;
}


void Timer::on()
{
	if(ons!=0)
	{
		static bool a;
		if(!a)cerr<<"Timings are unreliable due to nested timers."<<endl;
		a=true;
	}
//	assert(ons==0);
  ons++;
  n++;

  if(doTimingThisTime())
    {
      timerStartSec=time(0);
      timerStart=clock();
    }
}


void Timer::off()
{
  if(doTimingThisTime())
    {
      int t2=time(0);
      if((t2-timerStartSec)>(1000))
	total+= (t2-timerStartSec);
      else
	total+=((((double)clock()-(double)timerStart))*((1.0)/CLOCKS_PER_SEC));
    }
  ons--;
  assert(ons>=0);
}


void Timer::print(FILE *f)
{
  if(total!=0)
    fprintf(f,"%24s %8.2f s\n",s,(float)(total*skip));
}


void Timer::printList(FILE *f)
{
  fprintf(f,"Printing timer list:\n");
  Timer *l=timers;
  while(l)
    {
      l->print(f);
      l=l->next;
    }
}
