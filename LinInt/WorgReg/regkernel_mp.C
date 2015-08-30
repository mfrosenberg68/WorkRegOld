#include <math.h>
#include "mpmod.H"
#include "data_mp.H"

mp_real mpregkerne2(mp_real x_, mp_real x0_, mp_real b_)
{
 static mp_real pi = mppic ;
 mp_real _y = exp(-b_ * b_ * (x0_ -x_) *(x0_ - x_)) * cos(pi * b_ * b_ * (x0_ -x_ ));
 return _y;
}

mp_real mpkern(mp_real x_, mp_real x0_, mp_real b_, mp_real xi)
{
mp_real _y = mpe2(x_, xi) * mpregkerne2(x_, x0_, b_);
return _y;
}

mp_real mpkernpurereg(mp_real x_, mp_real x0_, mp_real b_)
{
mp_real _y = mpe2exakt(x_) * mpregkerne2(x_, x0_, b_);
return _y;
}

