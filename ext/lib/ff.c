// MolModExt implements a few number crunching routines for the molmod package in C.
// Copyright (C) 2007 - 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>
//
// This file is part of MolModExt.
//
// MolModExt is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// MolModExt is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>
//
// --


#include <math.h>
#include <stdlib.h>

inline void add_grad(int i, int j, double s, double *delta, double *gradient) {
  gradient[i*3  ] += s*(delta[0]);
  gradient[j*3  ] -= s*(delta[0]);
  gradient[i*3+1] += s*(delta[1]);
  gradient[j*3+1] -= s*(delta[1]);
  gradient[i*3+2] += s*(delta[2]);
  gradient[j*3+2] -= s*(delta[2]);
}

//takes locations of start and end point coordinate arrays (i and j) and location to write shortest vector array
inline void calc_delta(double *i, double *j, double *delta, double *unitcell, double *unitcell_reciproke, int *unitcell_active)
{
    double unit_cell[] = {1,0,0,0,2,0,0,0,1};
    double unit_cell_reciproke[] = {1,0,0,0,0.5,0,0,0,1};

    int p;
    for(p=0;p<3;p++)
    {
    if (*(unitcell_active+p) == 1)
        {//  minimal image convention - count distance between nearest 'images' of points

            int k,l;
            double tmp;
            double i_cart[3], j_cart[3], i_frac[3], j_frac[3], delta_frac[3];
            for(k=0;k<3;k++) // p coordinate (x,y or z) to fractional coordinates for both i and j
            {
                for(l=0;l<3;l++)
                {
                    i_cart[l] = *(i+l);
                    tmp = unit_cell[k+l]*i_cart[l];
                    i_frac[k] += tmp;
                    
                    j_cart[l] = *(j+l);
                    tmp = unit_cell[k+l]*j_cart[l];
                    j_frac[k] += tmp;
                }
            }

            //delta in frac coords
            for(k=0;k<3;k++)
            {
                delta_frac[k] = (j_frac[k] - i_frac[k])-round(j_frac[k] - i_frac[k]);
            }
            
            //delta back to cart for relevant coordinate
            double delta_cart_p = 0;
            for(k=0;k<3;k++) // p coordinate (x,y or z) to fractional coordinates for both i and j
            {
                tmp = unit_cell_reciproke[k+l]*delta_frac[k];
                delta_cart_p += tmp;
            }
            
            //write relevant delta coord to mem location
            *(delta+p) = delta_cart_p;
             //printf("Active periodic direction - delta with min. img. convention: %f\n", delta_frac);
        }
        else
        {
            *(delta+p) = *(i+p)-*(j+p);
            //printf("None-active direction - delta in standard way: %f\n", *(delta+p));
        }
    }

}

double ff_dm_quad(int n, double *cor, double *dm0, double *dmk, double amp, double *unitcell, double *unitcell_reciproke, int *unitcell_active, double *gradient) {
  int i,j;
  double d, d0, k, tmp, result;
  double delta[3];

  //printf("Unitcell active: %d %d %d \n",*(unitcell_active),*(unitcell_active+1),*(unitcell_active+2));
  result = 0.0;
  //printf("n=%i\n", n);
  for (i=0; i<n; i++) {
    for (j=0; j<i; j++) {
      d0 = dm0[i*n+j];
      k = dmk[i*n+j];
      //printf("i=%i  j=%i  d0=%i\n", i,j,d0);
      if (d0>0) {
        calc_delta(cor+3*i, cor+3*j, delta, unitcell, unitcell_reciproke, unitcell_active);

        tmp = delta[0];
        d = tmp*tmp;
        tmp = delta[1];
        d += tmp*tmp;
        tmp = delta[2];
        d += tmp*tmp;
        d = sqrt(d);

        //printf("d=%f    d0=%f\n", d, d0);
        //tmp = (d-d0*(radii[i]+radii[j]));
        //result += amp*tmp*tmp*radii[i]*radii[j]/d0;
        tmp = (d-d0);
        result += amp*k*tmp*tmp;
        if (gradient!=NULL) {
          //tmp = 2*amp*tmp/d0*radii[i]*radii[j]/d;
          tmp = 2*amp*k*tmp/d;
          add_grad(i, j, tmp, delta, gradient);
        }
        //result += tmp*tmp;
      }
    }
  }
  //printf("result=%f\n", result);
  //printf("delta=%f\n",delta[0]);
  return result;
}


double ff_dm_reci(int n, double *radii, double *cor, int *dm0, double amp, double *unitcell, double *unitcell_reciproke, int *unitcell_active, double *gradient) {
  int i,j,d0,r0;
  double d, tmp, result;
  double delta[3];

  result = 0.0;
  for (i=0; i<n; i++) {
    for (j=0; j<i; j++) {
      d0 = dm0[i*n+j];
      if (d0>1) {
        calc_delta(cor+3*i, cor+3*j, delta, unitcell, unitcell_reciproke, unitcell_active);

        tmp = delta[0];
        d = tmp*tmp;
        tmp = delta[1];
        d += tmp*tmp;
        tmp = delta[2];
        d += tmp*tmp;
        d = sqrt(d);
        r0 = radii[i]+radii[j];
        if (d < r0) {
            d /= r0;
            result += amp*(d-1)*(d-1)/d/d0;
            if (gradient!=NULL) {
              tmp = amp*(1-1/d/d)/r0/d/d0;
              add_grad(i, j, tmp, delta, gradient);
            }
        }
      }
    }
  }
  return result;
}


double ff_bond_quad(int m, int n, double *cor, int *pairs, double *lengths, double amp, double *unitcell, double *unitcell_reciproke, int *unitcell_active, double *gradient) {
  int b, i, j;
  double result, d, tmp;
  double delta[3];

  result = 0.0;
  //printf("m=%i\n", m);
  for (b=0; b<m; b++) {
    i = pairs[2*b  ];
    j = pairs[2*b+1];
        calc_delta(cor+3*i, cor+3*j, delta, unitcell, unitcell_reciproke, unitcell_active);

        tmp = delta[0];
        d = tmp*tmp;
        tmp = delta[1];
        d += tmp*tmp;
        tmp = delta[2];
        d += tmp*tmp;
        d = sqrt(d);

    tmp = d-lengths[b];
    result += amp*tmp*tmp;
    if (gradient!=NULL) {
      tmp = 2*amp*tmp/d;
      add_grad(i, j, tmp, delta, gradient);
    }
    //printf("result=%f\n", result);
  }
  return result;
}

double ff_bond_hyper(int m, int n, double *cor, int *pairs, double *lengths, double scale, double amp, double *unitcell, double *unitcell_reciproke, int *unitcell_active, double *gradient) {
  int b, i, j;
  double result, d, tmp;
  double delta[3];

  result = 0.0;
  for (b=0; b<m; b++) {
    i = pairs[2*b  ];
    j = pairs[2*b+1];
    calc_delta(cor+3*i, cor+3*j, delta, unitcell, unitcell_reciproke, unitcell_active);

    tmp = delta[0];
    d = tmp*tmp;
    tmp = delta[1];
    d += tmp*tmp;
    tmp = delta[2];
    d += tmp*tmp;
    d = sqrt(d);

    tmp = d-lengths[b];
    result += amp*(cosh(scale*tmp)-1);
    if (gradient!=NULL) {
      tmp = amp*scale*sinh(scale*tmp)/d;
      add_grad(i, j, tmp, delta, gradient);
    }
  }
  return result;
}


