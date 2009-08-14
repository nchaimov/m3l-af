/* randist/pareto.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
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

#include <config.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* The Pareto distribution has the form,

   p(x) dx = (a/b) / (x/b)^(a+1) dx     for x >= b

 */

double
gsl_ran_pareto (const gsl_rng * r, double a, const double b)
{
  double x = gsl_rng_uniform_pos (r);

  double z = pow (x, -1 / a);

  return b * z;
}

double
gsl_ran_pareto_pdf (const double x, const double a, const double b)
{
  if (x >= b)
    {
      double p = (a/b) / pow (x/b, a + 1);
      return p;
    }
  else
    {
      return 0;
    }
}
