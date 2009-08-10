/* specfunc/test_sf.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Author:  G. Jungman */

#ifndef TEST_SF_H
#define TEST_SF_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_sf_result.h>

double test_sf_frac_diff(double x1, double x2);
int test_sf_check_result(char * message_buff, gsl_sf_result r, double val, double tol);
int test_sf_check_val(char * message_buff, double rval, double val, double tol);
int test_sf_check_return(char * message_buff, int val_return, int expected_return);
int test_sf_check_result_relax(char * message_buff, gsl_sf_result r, double val, double tol);

#define TEST_TOL0  (2.0*GSL_DBL_EPSILON)
#define TEST_TOL1  (16.0*GSL_DBL_EPSILON)
#define TEST_TOL2  (256.0*GSL_DBL_EPSILON)
#define TEST_TOL3  (2048.0*GSL_DBL_EPSILON)
#define TEST_TOL4  (16384.0*GSL_DBL_EPSILON)
#define TEST_TOL5  (131072.0*GSL_DBL_EPSILON)
#define TEST_TOL6  (1048576.0*GSL_DBL_EPSILON)
#define TEST_SQRT_TOL0 (2.0*GSL_SQRT_DBL_EPSILON)
#define TEST_SNGL  (1.0e-06)

#define TEST_SF_INCONS  1
#define TEST_SF_ERRNEG  2
#define TEST_SF_TOLBAD  4
#define TEST_SF_RETBAD  8

int test_sf (gsl_sf_result r, double val_in, double tol, int status, int expect_return, const char * desc);
int test_sf_val (double val, double val_in, double tol, const char * desc);
int test_sf_rlx (gsl_sf_result r, double val_in, double tol, int status, int expect_return, const char * desc);
int test_sf_2 (gsl_sf_result r1, double val1, double tol1, gsl_sf_result r2, double val2, double tol2, int status, int expect_return, const char * desc);
int test_sf_sgn (gsl_sf_result r, double sgn, double val_in, double tol, double expect_sgn, int status, int expect_return, const char * desc);

#define TEST_SF(stat, func, args, val_in, tol, expect_return) { int status = func args; stat += test_sf(r, val_in, tol, status, expect_return, #func #args); }

#define TEST_SF_VAL(stat, func, args, val_in, tol) { double val = func args; stat += test_sf_val(val, val_in, tol, #func #args); }

#define TEST_SF_RLX(stat, func, args, val_in, tol, expect_return) { int status = func args; stat += test_sf_rlx(r, val_in, tol, status, expect_return, #func #args); }

#define TEST_SF_2(stat, func, args, val1, tol1, val2, tol2, expect_return) { int status = func args; stat += test_sf_2(r1, val1, tol1, r2, val2, tol2, status, expect_return, #func #args); }

#define TEST_SF_SGN(stat, func, args, val_in, tol, expect_sgn, expect_return) { int status = func args; stat += test_sf_sgn(r, sgn, val_in, tol, expect_sgn, status, expect_return, #func #args); }

int test_airy(void);
int test_bessel(void);
int test_coulomb(void);
int test_dilog(void);
int test_gamma(void);
int test_hyperg(void);
int test_legendre(void);


#endif /* !TEST_SF_H */
