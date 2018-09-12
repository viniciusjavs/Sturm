#include "sturm.h"

int main()
{
    Sturm sturm;
    sturm.set_parameters();
    sturm.show_sturm_sequence(sturm.get_sturm_sequence());
    sturm.show_roots(sturm.get_real_roots());
}
