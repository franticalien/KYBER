#ifndef INVERSE_H
#define INVERSE_H

#include <stdint.h>
#include <stdio.h>
#include "params.h"
#include "polyvec.h"
#include "poly.h"

#define mat_inv_2 KYBER_NAMESPACE(mat_inv_2)
void mat_inv_2(polyvec *b, polyvec *a);

#endif


