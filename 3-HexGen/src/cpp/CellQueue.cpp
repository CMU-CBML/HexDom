//------------------------------------------------------------------------
//
// cellQueue.C - queue of cell identifiers.  The circular queue dyanmically
//               resizes itself when full.  Elements in the queue are of
//               indicies of i,j,k (specialized for 3d structured grids)
//
// Copyright (c) 1997 Dan Schikore
//------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <malloc.h>
#include <string.h>
#include "cellQueue.h"
