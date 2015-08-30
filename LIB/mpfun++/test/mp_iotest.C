#include "mpmod.H"

int main()

{
    mp_init();
    ifstream from("mp_iotest.in");
    if (! from)
      cerr << "Couldn't open input file\n";
    ofstream to("mp_iotest.out");
    if (! to)
      cerr << "Couldn't open output file\n";
    mp_real r;
    while (from >> r) to << r;
    return 0;
}
