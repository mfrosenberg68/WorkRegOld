/*
%  Dortmund C++ Class and Template Library 
%  Copyright (C) 1994 Wolfgang Wenzel and others

%  The code in this file was derived from routines published in 
%  Numerical Recipes, it its the responsibility of the user to insure that
%  he is authorized to use them. Please read the appropriate section entitled
%  LICENSE in the documentation. 

%  This library is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%  Library General Public License for more details.

%  You should have received a copy of the GNU Library General Public
%  License along with this library; if not, write to the Free
%  Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%
% $Id: xtime.c,v 1.2 1995/10/06 15:05:46 wenzel Exp $
%
*/

#include <iostream.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "xtime.h"

Timer* timer;

Timer::Timer(const int mx,const int mxdep)
{
  max      = mx;
  maxdepth = mxdep;  

  msec  = new  long[max];
  assert(msec !=0);
  count = new  long[max];
  assert(count !=0);
  name  = new char*[max];
  assert(name != 0);
  stime = new  long[maxdepth];
  assert(stime !=0);
  code = new   int[maxdepth];
  assert(code !=0);
  
  for(int i=0; i<max;i++)
  {
    msec[i]  = 0;
    count[i] = 0;
    name[i]  = strdup("FUNCTION");
  }
  level = 0;  
}

Timer::~Timer() 
{
  delete msec;
  delete count;
  delete stime;  
  for(int i=0; i < max; i++)
    delete name[i];
  delete name;
}

void Timer::print()
{
  char buffer[100];
  
  cout << " +++ TIMING INFORMATION: " << endl;
  cout << endl << "               NAME   NO      COUNT     TIME       TIME/CALL " << endl;
  for(int i=0;i<max;i++)
    if (count[i] > 0) 
    {
      double ms = get_time(i);
      sprintf(buffer,"    %15s %4i %10i %10.3f    %10.3f \n",
	      name[i],i,count[i],ms/1E3,1.0*ms/count[i]);
      cout <<  buffer;
    }
  cout << endl;
  
}





