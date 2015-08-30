#ifdef xlC
#define LANG   "FORTRAN"
#else
#define LANG   "C"
#endif

#ifdef ibmrs6000
#define SIGN   sign
#define MPADD  mpadd
#define MPANG  mpang
#define MPCADD mpcadd
#define MPCDIV mpcdiv
#define MPCEQ  mpceq
#define MPCMUL mpcmul
#define MPCPR  mpcpr
#define MPCPWR mpcpwr
#define MPCSQR mpcsqr
#define MPCSSH mpcssh
#define MPCSSN mpcssn
#define MPCSUB mpcsub
#define MPDIV  mpdiv
#define MPDMC  mpdmc
#define MPEQ   mpeq
#define MPEXP  mpexp
#define MPINFR mpinfr
#define MPINPC mpinpc
#define MPLOG  mplog
#define MPMDC  mpmdc
#define MPMMPC mpmmpc
#define MPMUL  mpmul
#define MPMULD mpmuld
#define MPNINT mpnint
#define MPNPWR mpnpwr
#define MPNRT  mpnrt
#define MPOUT  mpout
#define MPOUTS mpouts
#define MPPI   mppi
#define MPRAND mprand
#define MPSETP mpsetp
#define MPSQRT mpsqrt
#define MPSUB  mpsub
#else
#define SIGN   r_sign
#define MPADD  mpadd_
#define MPANG  mpang_
#define MPCADD mpcadd_
#define MPCDIV mpcdiv_
#define MPCEQ  mpceq_
#define MPCMUL mpcmul_
#define MPCPR  mpcpr_
#define MPCPWR mpcpwr_
#define MPCSQR mpcsqr_
#define MPCSSH mpcssh_
#define MPCSSN mpcssn_
#define MPCSUB mpcsub_
#define MPDIV  mpdiv_
#define MPDMC  mpdmc_
#define MPEQ   mpeq_
#define MPEXP  mpexp_
#define MPINFR mpinfr_
#define MPINPC mpinpc_
#define MPLOG  mplog_
#define MPMDC  mpmdc_
#define MPMMPC mpmmpc_
#define MPMUL  mpmul_
#define MPMULD mpmuld_
#define MPNINT mpnint_
#define MPNPWR mpnpwr_
#define MPNRT  mpnrt_
#define MPOUT  mpout_
#define MPOUTS mpouts_
#define MPPI   mppi_
#define MPRAND mprand_
#define MPSETP mpsetp_
#define MPSQRT mpsqrt_
#define MPSUB  mpsub_
#endif

#ifdef doubleonly
#define FHYPOT hypot
#define ATAN2F atan2
#define COSF   cos
#define SINF   sin
#define FEXP   exp
#define POWF   pow
#define SQRTF  sqrt
#else
#define FHYPOT fhypot
#define ATAN2F atan2f
#define COSF   cosf
#define SINF   sinf
#define FEXP   fexp
#define POWF   powf
#define SQRTF  sqrtf
#endif
