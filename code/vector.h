#ifndef VECTOR_H
#define VECTOR_H

#include "main.h"
#include "vector.cpp"

/**
 * Calculate the norm of a 3d vector.
 */
double v3_norm(const listdouble& v3);

/**
 * Calculate the cross product of two 3d vectors.
 */
listdouble v3_cross(const listdouble& v1, const listdouble& v2);

/**
 * Calculate the scalar product of two 3d vectors.
 */
double v3_scalar(const listdouble& v1, const listdouble& v2);

#endif