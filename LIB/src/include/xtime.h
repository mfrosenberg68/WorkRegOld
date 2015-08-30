#ifndef TIMER_DEF
#define TIMER_DEF
#include <time.h>
#include <assert.h>
#include <string.h>

#ifdef GNU
const long timer_scale = 1;
#else
const long timer_scale = 1000;
#endif

const long timer_fix   = 4294967;

class Timer
{
  int max;
  char         **name;
  long *msec;
  long *count;  
  long *stime;
  int  *code;
  int level;
  int maxdepth;
  
public:
  Timer(const int mx=100,const int mxdep = 20);
 ~Timer();

long get_time(int i) const 
{ return (long) (1000.0*msec[i]/CLOCKS_PER_SEC*timer_scale); } 
void set_name(int i,char* newname)
{
  assert(i >=0 && i < max);
  delete name[i];
  name[i] = strdup(newname);
}

void start(const int timer_code )  
{ 
  code [++level] = timer_code;
  stime[  level] = clock() / timer_scale;
  if (level >= maxdepth) cerr << "Timer Overflow." << endl;
}
void stop(const int timer_code )   
{
  if (timer_code != code[level])
  {
    cout << "Error in Timer: Trying to Stop: " << timer_code 
      << " expected: " << code[level] << endl;
  }
  
  long etime = (long) clock()/ timer_scale;
  if (etime < stime[level])
  {
    cout << endl << "TIMER OVERFLOW level: " << level
         << " START: " << stime[level] << " END: " << etime << " NEW: ";
    for(int i =0 ; i < level; i++)
    {
      msec[i]  += etime + timer_fix - stime[i];
      stime[i]  = etime;
    }
    etime += timer_fix;
    cout << etime << " DIFF: " << etime - stime[level] <<  endl;
  }
  msec [timer_code] += etime - stime[level];  
  count[timer_code]++;
  if (--level<0) 
    cerr << "Timer Underflow" << endl; 
}
void print();
};

extern Timer* timer;

#endif
