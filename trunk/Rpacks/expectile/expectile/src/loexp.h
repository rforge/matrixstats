/*
 * loexp - find expectile curves of local polynomial regression 
 *
 *  Copyright (C) 1999-2008, Pratyaksha J. Wirapati
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available via WWW at
 *  http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
 *  writing to the Free Software Foundation, Inc., 51 Franklin Street
 *  Fifth Floor, Boston, MA 02110-1301  USA.

*/

// 
// loexp() fits local expectiles.
//
// The input is series of (x0,y0),(x1,y1),...  However, x_i - x_{i-1}
// is assumed constant, and the x values are not considered. Each observation
// can be given non-negative weight. This routine find 
//   E(y) = beta0 + beta1 x + beta2 x^2
// at fitted locally each point, using Gaussian kernel to weigh nearby
// points. The kernel width is determined by sigma (defined in the number
// of data points).
//
// beta's for each data point i were found to minimize
//
//   \sum_{j=1}^n w_j G_{ij} R_j A_j (y_j-E[y|x_j,beta_i])^2
// 
// which is local (G is gaussian kernel), robust (R is Tukey's biweight 
// determined by the median of residuals) and asymmetric (A is the expectile
// weight: it is 1-alpha if the residual is negative and alpha if positive).
// 
// The robustness and asymmetry weights are pooled across all data points
// (they are not local). However, the robustness weights are separate fo
// positive and negative residuals.
// 
// The polynomial dummy x's is arbitrary and thus beta's are not
// meaningful and not returned. Ey is invariant to the choice of dummies.
//
// return:
//  0  converge
//  1  iteration limit exceeded
//  2  non-positive definite X'WX encoutered when solving local regression
//  3  invalid order of polynomials
//  
int
loexp (
    int n,                // data length
    const double *y,      // data vector
    const double *w,      // sample weights (assume w[i]=1 if w == NULL)
    int polynomial_order, // 0, 1 or 2
    double sigma,         // gaussian kernel SD (in # of data points)
    double alpha,         // asymmetry weight [0..1]
    double biweight,      // Tukey's biweight tuning

    double tolerance,     // convergence tolerance
    int iter_max,         // maximum iteration
    
    // output (arrays should be allocated by caller)
    double *Ey,           // fitted curve
    double *v_            // final weights (biweight * asymmetry), can be NULL
    )
  ;
