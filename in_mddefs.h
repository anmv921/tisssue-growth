#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <stdio.h> 
#include <algorithm>

typedef long double real;

#include "in_vdefs.h"
#include "in_proto.h"

#define DO_CELLS for (n = 0; n < Ncells; n ++)

#define VWrap(v, t) \
   if (v.t >= 0.5 * region.t)      v.t -= region.t; \
   else if (v.t < -0.5 * region.t) v.t += region.t

#define VWrapAll(v) \
   {VWrap (v, x); \
   VWrap (v, y);}
