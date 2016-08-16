#ifndef TIMER_H_INCLUDED
#define TIMER_H_INCLUDED

#include <stdio.h>

class Timer
{
  static class Timer *timers;
  class Timer *next;
  const char *s;
  int skip;
  int ons;
  int n;
  double total;
  int timerStart;
  int timerStartSec;
  bool doTimingThisTime();
 public:
  Timer(const char* s, int skip=10);
  void on();
  void off();
  void print(FILE *f);
  static void printList(FILE *f=stderr);
};

#if 1

class TimerScope
{
  Timer *timer;
 public:
  TimerScope(Timer *timer){this->timer=timer;timer->on();}
  ~TimerScope(){timer->off();}
};

#else

class TimerScope
{
 public:
  TimerScope(Timer *timer){}
  ~TimerScope(){}
};

#endif

#endif
