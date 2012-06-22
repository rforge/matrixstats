/*
sfit.c 
*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>


// robustness weights

static inline double
w_huber ( double absres, double par )
{
  return absres < par ? 1 : par/absres ;
}

static inline double
w_biweight ( double absres, double par )
{
  if( absres > par ) return 0;
  double u = absres/par;
  u = 1-u*u;
  u *= u;
  return u;
}

static inline void
update_moment3 ( double n, double M1, double M2, double M3, double x,
  double *n_, double *M1_, double *M2_, double *M3_ )
{
  if(!isfinite(x)) return;
  *n_ = ++n;
  double d = (x-M1)/n;
  *M1_ = M1 + d;
  *M2_ = M2 + n*(n-1)*d*d;
  *M3_ = M3 + n*(n-1)*(n-2)*d*d*d - 3*M2*d;
}


// XXX calling wquantile(0,...) results in access to x[-1]

// Find a quantile of weighted data
// 
// input:
// x[0],...,x[n-1] and corresponding non-negative weights w[0],...,w[n-1]
// The data is not sorted (otherwise the problem is trivial)
// 
// return x[k], such that sum_i^k w[i] <= p < sum_i^{k+1} w[i]
// The indices are such that the cumulative sums correspond to those
// of sorted data (although only partial quicksort is done here). 
// Both x and w will be scrambled by this routine.
// 
// example: p = sum_w/2 is the weighted median
//
// notes:
// - the weights don't have to sum to one; but if it does, p is the
//   cumulative probability
// - the function x(p), where p is the cumulative sum of weights,
//   is assumed to be a step function defined on intervals
//      [-inf,p_x0),[p_x0,p_x1),[p_x1,p_x2),...[p_x{n-1},+inf]
//   Thus, the larger x is returned if p is exactly on an interval boundary
// 
static double
wquantile ( 
    int n, double *x, // data
    double *w,        // weights corresponding to x, w[i] >= 0
    double sum_w,     // sum w[i]
    double p          // requested cumulative sum
    )
{
  double t, v;
  int a = 0, b = n-1;
  double sa = 0, sb = sum_w;
  while( b > a )
    {
    // put the midpoint of three at b
    if( x[b] < x[a] )
      {
      t = x[a]; x[a] = x[b]; x[b] = t; 
      v = w[a]; w[a] = w[b]; w[b] = v; 
      }
    int c = (a+b)/2;
    if( x[c] < x[b] )
      {
      t = x[c]; x[c] = x[b];
      v = w[c]; w[c] = w[b];
      if( t < x[a] )
        {
        x[b] = x[a]; x[a] = t;
        w[b] = w[a]; w[a] = v;
        }
      else
        { x[b] = t; w[b] = v; }
      }
    
    int i = a-1, j = b;
    double sl = sa, sr = sb;
    for(;;)
      {
      while( x[++i] < x[b] )
        sl += w[i];
      while( x[b] < x[--j] )
        sr -= w[j];
      if( i >= j ) break;
      t = x[i]; x[i] = x[j]; x[j] = t;
      v = w[i]; w[i] = w[j]; w[j] = v;
      sl += w[i]; sr -= w[j];
      }
    t = x[i]; x[i] = x[b]; x[b] = t;
    v = w[i]; w[i] = w[b]; w[b] = v;
    
    if( p < sl )
      { b = i - 1; sb = sl; }
    else if( p >= sr )
      { a = i + 1; sa = sr; }
    else
      return x[i];
    }
  return x[b];
}

#if 1

/* Choose Householder vector 
 *
 * x[1 : n-1] is overwritten by the Householder vector u (which is scaled
 * such that the first element is 1)
 * x[0] is overwritten by the norm of x (which is also the result of H_u x)
 * return: tau = 2/(u' u)
 * 
 * ref: Golub & van Loan (1996) "Matrix Computations" 2nd ed., chapter 5
 */

static inline double
hoho_vector ( int n, double *x )
{
  double ss = 0;
  for(int i = 1; i < n; i++ )
    ss += x[i]*x[i];
  if( ss < DBL_MIN )
    return 0;
    
  double a, normx = sqrt( x[0]*x[0] + ss );
  if( x[0] <= 0 )
    a = x[0] - normx;
  else
    a = - ss /( x[0] + normx );
  double tau = 2*a*a/(a*a+ss);
  
  x[0] = normx ;
  for(int i = 1; i < n; i++ )
    x[i] /= a;
  return tau;
}

/*
 * apply Householder transformation (I - tau u u') x 
 * x is overwritten.
 */

static inline void
hoho_reflect ( int n, double tau, const double *u, double *x )
{
  if( tau < DBL_MIN ) return;
  double b = x[0];
  for( int j = 1; j < n; j++ )
    b += u[j] * x[j];
  b *= tau;
  x[0] -= b;
  for( int j = 1; j < n; j++ )
    x[j] -= b * u[j];
}

// as above, but apply to a submatrix
static inline void
hoho_reflect2 (
  int m,              /* number of vectors */
  int n,              /* length of each vector */
  double tau,         /* 2/(u' u), as produced by hoho_vector */
  const double *u,    /* householder vector, as produced by hoho_vector */
  double *x, int dx   /* dx is displacement to the next vector */
  )
{
  if( tau < DBL_MIN ) return;
  for(int i = 0; i < m; i++, x += dx )
    {
    double b = x[0];
    for( int j = 1; j < n; j++ )
      b += u[j] * x[j];
    b *= tau;
    x[0] -= b;
    for( int j = 1; j < n; j++ )
      x[j] -= b * u[j];
    }
}

/* solve by back substitution */
static inline void
backsub (
  int m,                      // number of terms 
  const double *R, int dR,    // R as produced by hoho_QR
  double tolerance,           // smallest 
  double *b )                 // right-hand-side (overwritten)
{
  R += (m-1)*dR;
  for( int i = m-1; i >= 0; i-- )
    {
    double sum = b[i]; const double *r = R + i;
    for( int j = m - 1; j > i; j--, r -= dR )
      sum -= b[j] * r[0];
    if( fabs(r[0]) > tolerance )
      b[i] = sum / r[0];
    else
      b[i] = 0;
    }
}
#endif

//
// Find initial simplex.
// 
void
initial_simplex (
    int n, int p, const double *Y,  // data matrix
    int m,            // number of vertices
    double shrink,    // shrinkage towards the mean of vertices, [0..1)
    double *X         // initial simplex   
    )   
{
  // permutation matrix for selecting vertices from data points
  int* P = (int*)malloc(sizeof(int)*n);
  for(int j = 0; j < n; j++ )
    P[j] = j;
  
  // find the mean
  double *workspace = (double*)malloc(sizeof(double)* (p + p + p*n + n));
  double *mu = workspace;
  double *mu0 = mu + p;
  double *Z = mu0 + p;
  double *norm = Z + p*n;

  for(int k = 0; k < p; k++ )
    mu[k] = 0;
  const double *Yj = Y;
  for(int j = 0; j < n; j++, Yj += p )
    for(int k = 0; k < p; k++ )
      mu[k] += Yj[k];
  for(int k = 0; k < p; k++ )
    mu[k] /= n;
  
  // find the furthest point from the mean
  int maxj = -1; double maxnorm = 0;
  Yj = Y;
  for(int j = 0; j < n; j++, Yj += p )
    {
    double normj = 0;
    for(int k = 0; k < p; k++ )
      normj += (Yj[k]-mu[k])*(Yj[k]-mu[k]);
    if( normj > maxnorm )
      { maxnorm = normj; maxj = j; }
    }
   
  // use as the first vertex
  P[maxj] = 0; P[0] = maxj;

  // copy the centered, Y[j]-X[0], points. compute and save their norms
  // select the second vertex as the one with largest residual norm
  // when projected to mean-X[0]
  // 
  maxj = -1; maxnorm = 0;
  const double *X0 = Y + P[0]*p; 
  double sum = 0;
  for(int k = 0; k < p; k++ )
    {
    mu0[k] = mu[k] - X0[k];
    sum += mu0[k]*mu0[k];
    }
  for(int k = 0; k < p; k++ )
    mu0[k] /= sqrt(sum);

  double *Zj = Z;
  Yj = Y;
  for(int j = 0; j < n; j++, Zj += p, Yj += p )
    {
    norm[j] = 0;
    sum = 0;
    for(int k = 0; k < p; k++ )
      {
      Zj[k] = Yj[k] - X0[k];
      norm[j] += Zj[k]*Zj[k];
      sum += Zj[k]*mu0[k];
      }
    double resnorm = norm[j]-sum*sum;
    if( resnorm > maxnorm )
      { maxnorm = resnorm; maxj = j; }
    }

  // QR decomposition with pivoting (depends on maxnorm & maxj above)
  for(int i = 0; i < m-1; i++ )
    {
    int t = P[i+1]; P[i+1] = maxj; P[maxj] = t;

    double *zi = Z + maxj*p;
    double tau = hoho_vector( p-i, zi + i );
    maxj = -1; maxnorm = 0;
    for(int a = i+1; a < n-1; a++ )
      {
      int j = P[a+1];
      Zj = Z + j*p;
      hoho_reflect( p-i, tau, zi + i, Zj + i );
      double v = Zj[i];
      norm[j] -= v*v;
      if(norm[j] > maxnorm )
        { maxnorm = norm[j]; maxj = j; }
      }
    }
  
  for(int i = 0; i < m; i++ ) // copy and compute the mean
    {
    double *xi = X + i*p;
    const double *yi = Y + P[i]*p;
    for(int k = 0; k < p; k++ )
      xi[k] = yi[k];
    }

  // find the apex
  double maxskew = -INFINITY; int imaxskew = -1;
  for(int i = 0; i < m; i++ )
    {
    double *xi = X + i*p;
    for(int k = 0; k < p; k++ )
      Z[k] = mu[k]-xi[k];
    double N = 0, M1 = 0, M2 = 0, M3 = 0;
    Yj = Y;
    for(int j = 0; j < n; j++, Yj += p )
      {
      double dot = 0;
      for(int k = 0; k < p; k++ )
        dot += Yj[k]*Z[k];
      update_moment3 ( N, M1, M2, M3, dot, &N, &M1, &M2, &M3 );
      }
    M2 /= N-1;
    double skewness = M3/N / sqrt(M2*M2*M2);
    if( skewness > maxskew )
      { maxskew = skewness; imaxskew = i; }
    }
    
  // put the apex in X[0]
  for(int k = 0; k < p; k++ )
    {
    double t = X[k]; X[k] = X[imaxskew*p+k]; X[imaxskew*p+k] = t;
    }

  if(shrink > 0 && shrink < 1)
    for(int i = 0; i < m; i++ )
      {
      double *xi = X + i*p;
      for(int k = 0; k < p; k++ )
        xi[k] = (1-shrink)*xi[k] + shrink*mu[k];
      }

  free(workspace); free(P);
}



int
sfit0 (
  int n,                // data length
  int p,                // data dimensions
  const double *Y,      // data matrix, n x p
  const double *w,      // observation weight, n
  int m,                // simplex order
  int autoinit,         // 1: call initial_simplex, 0: use caller's X
  double *X,            // initial simplex (overwritten by the final), m x p
  double *wX0,          // prior weights
  double *X0,           // prior simplex/cone (depending on autoinit)
  double *Beta,         // affine combination coefficients, n x m

  double lambda,        // vertex assigment parameter
  double alpha,         // asymmetry weight
  int M_family,   // 0: least-squares, 1: Huber, 2: Biweight
  double gamma,         // robustness scale factor
                        // if < 0, use 1.345 for Huber and 4.685 for biweight
  int fitcone,          // if 1, treat first vertex as the apex (use separate
                        // robust scale estimate for the opposite facet).
                        // Otherwise, treat as simplex and pool scale est.
  
  int verbose,          // dump iteration steps to stderr
  double tolerance,     // convergence tolerance, ||\Delta X||F/||X||_F
  int iter_max,         // maximum iterations
  double Rtolerance     // smallest acceptable abs value of the diagonal entry
                        // of R still used for solving. 
  )
{
  if( m == 0 || n == 0 || p == 0 ) return 0;

  int prior_simplex = 0;
  if( wX0 && X0 )
    {
    double sum_wX0 = 0;
    for(int i = 0; i < m; i++ )
      {
      if( wX0[i] < 0 ) wX0[i] = 0;
      sum_wX0 += wX0[i];
      }
    if( sum_wX0 > 0 )
      {
      prior_simplex = 1;
      if( fitcone == 1 ) // add the apex coordinate to the others
        {
        double *X0i = X0+p;
        for(int i = 1; i < m; i++, X0i += p )
          for(int j = 0; j < p; j++ )
            X0i[j] += X0[j];
        }
      }
    }

  if(M_family > 0 && gamma < 0)
    gamma = (M_family == 1 ? 1.345 : 4.685 );
  fitcone = (fitcone == 0 ? 0: 1);
  
  double *workspace = (double*)malloc(
      sizeof(double)*((p+1)*(1+m+n) + 2*m*m + m + 4*n) );
  
  double *ej = workspace;        // residual vector (temporary)
  double *Xc = ej + p + 1;      // centered simplex (use X[0] as origin), QR
  double *Yc = Xc + m*(p+1);    // centered data (and Q'Yc and (I-Q')Yc
  double *D = Yc + n*(p+1);     // matrix of delta_{ik}
  double *sum_Wd = D + m*m;     // sums of weight for delta (vertex nudge)
  double *sum_We = sum_Wd + m*m; // sums of weight for residual twist
  double *negbeta = sum_We + m;    // for finding MAR
  double *wnegbeta = negbeta + n;  // weight for finding MAR
  double *e2 = wnegbeta + n;       // residual norm
  double *we2 = e2 + n;            // weight of residual
  
  if( autoinit )
    initial_simplex( n, p, Y, m, 0.1, X );
    
  for(int r = 0; r < iter_max; r++ )
    {
    
    // copy and center the current simplex 
    double *Xci = Xc + p + 1, *Xi = X + p;
    for(int i = 1; i < m; i++, Xci += p+1, Xi += p )
      for(int k = 0; k < p; k++ )
        Xci[k+1] = Xi[k] - X[k];

    // QR decompose Xc 
    // 
    // If a vertex collapsed (become absorbed into the opposite facet),
    // then the diagonal entry of R will be very small. This is handled
    // when solving for beta. beta_i is simply zeroed for all j's.
    // The consequence is that the vertex will not be updated.
    // This may be transient, as the simplex may be full-ranked again
    // in subsequent iterations. However, if this is due to intrinsic
    // rank deficiency in the data, the vertex will be
    // ignored until convergence. \sum_j w_j G(beta_ij,\lambda) can be used
    // as an indicator of support for vertex i.
    // 
    Xci = Xc + p + 1;
    for(int i = 1; i < m; i++, Xci += p+1 )
      {
      Xci[0] = hoho_vector( p+1-i, Xci+i );
      hoho_reflect2 ( m-1-i, p+1-i, Xci[0], Xci+i, Xci+i+p+1, p+1 );
      }
    
    int nneg = 0, nbase = n;
    double sumwneg = 0, sumwbase = 0, sumwe2 = 0;
    double *Ycj = Yc;
    const double *Yj = Y;
    for(int j = 0; j < n; j++, Ycj += p+1, Yj += p )
      {
      // copy original data
      for(int k = 0; k < p; k++ )
        Ycj[k+1] = Yj[k] - X[k];
      
      // Q' Yc[j]
      Xci = Xc + p + 1;
      for(int i = 1; i < m; i++, Xci += p+1 )
        hoho_reflect ( p+1-i, Xci[0], Xci+i, Ycj+i );
      
      // Solve R beta = Q'Yc[j], overwrite the R part of Yc
      backsub ( m-1, Xc+p+1+1, p+1, Rtolerance, Ycj+1 );

     // satisfy affine condition
      Ycj[0] = 1.0;
      for(int i = 1; i < m; i++ )
        Ycj[0] -= Ycj[i];

      // copy the most negative beta for MAR estimation
      // XXX option to treat beta_1j differently? (first vertex is an apex)
      double v = 0; int iv = -1;
      for(int i = 0; i < m; i++ )
        if( Ycj[i] < 0 && Ycj[i] < v )
          { v = Ycj[i]; iv = i; }

      if( v < 0 )
        {
        if( fitcone && iv == 0)
          {
          nbase--;
          negbeta[nbase] = -v;
          wnegbeta[nbase] = (w ? w[j]:1);
          sumwbase += wnegbeta[nbase];
          }
        else
          {
          negbeta[nneg] = -v;
          wnegbeta[nneg] = (w ? w[j]: 1);
          sumwneg += wnegbeta[nneg];
          nneg++;
          }
        }
      
      // compute residual norm, overwrite the storage
      if( p > m-1 )
        {
        double e = Ycj[m]*Ycj[m];
        for(int k = m+1; k <= p; k++ )
          e += Ycj[k]*Ycj[k];
        Ycj[m] = sqrt(e);
        e2[j] = Ycj[m];
        we2[j] = (w ? w[j]: 1);
        sumwe2 += we2[j];
        }
      }

    // get MAR and parameter for biweight robustness
    double bwbeta = gamma 
            * wquantile(nneg, negbeta, wnegbeta, sumwneg, sumwneg/2)/0.6745;
    double bwbase = 0;
    if( fitcone )
      bwbase = gamma * wquantile(n-nbase, negbeta+nbase, wnegbeta+nbase,
                                 sumwbase, sumwbase/2)/0.6745;

    double bwe2 = gamma 
            * wquantile(n, e2, we2, sumwe2, sumwe2/2)/0.6745;
    
       
    //
    // accumulate vertex nudge and residual twist
    //
    for(int i = 0; i < m*m; i++ )
      D[i] = sum_Wd[i] = 0;
    for(int i = 0; i < m; i++ )
      sum_We[i] = 0;
    for(int i = 0; i < m*p; i++ ) // overwrite by new X, use p as offset!
      Xc[i] = 0;
    
    Ycj = Yc;
    for(int j = 0; j < n; j++, Ycj += p+1 )
      {
      double *beta = Ycj;
      
      // vertex proximity weight
      double G[m], sum = 0;
      for(int i = 0; i < m; i++ )
        {
        G[i] = (beta[i] >= 1 ? 1: (beta[i] <= 0 ? 0: pow(beta[i],lambda)));
        sum += G[i];
        }

      for(int i = 0; i < m; i++ )
        G[i] /= sum;

      for(int k = 0; k < m; k++ )
        {
        double rpar = (fitcone && k == 0 ? bwbase: bwbeta );
        double u = 1;
        if( beta[k] < 0 && rpar > 0 )
          switch( M_family )
            {
            case 0:
              break;
            case 1:
              u = w_huber( -beta[k], rpar);
              break;
            case 2:
            default:
              u = w_biweight( -beta[k], rpar);
              break;
            }
          
        u *= (beta[k] <= 0 ? 1 - alpha : alpha ); // asymmetry
        u *= (w ? w[j]: 1);                       // observation
        for(int i = 0; i < m; i++ )
          {
          if( i == k ) continue;
          D[i*m+k] += u * G[i] * beta[k];
          sum_Wd[i*m+k] += u * G[i];
          }
        }
        
      if( p > m-1 ) // accumulate residual
        {
        Yj = Y + j*p;
        for(int k = 0; k < p; k++ )
          ej[k] = Yj[k];
        Xi = X;
        for(int i = 0; i < m; i++, Xi += p )
          for(int k = 0; k < p; k++ )
            ej[k] -= beta[i] * Xi[k];
        
        double u = 1;
        if( bwe2 > 0 )
          switch( M_family )
            {
            case 0:
              break;
            case 1:
              u = w_huber( beta[m], bwe2);
              break;
            case 2:
            default:
              u = w_biweight( beta[m], bwe2);
              break;
            }
        u *= (w? w[j] : 1); // observation
        
        Xci = Xc;
        for(int i = 0; i < m; i++, Xci += p )
          {
          double v = u * G[i];
          sum_We[i] += v;
          for(int k = 0; k < p; k++ )
            Xci[k] += v * ej[k];
          }
        }
      }

    // solve vertex nudge
    for(int i = 0; i < m; i++ )
      {
      double sum = 0;
      for(int k = 0; k < m; k++ )
        {
        if( k == i ) continue;
        if( sum_Wd[i*m+k] > 0 )
          D[i*m+k] /= sum_Wd[i*m+k];
        sum += D[i*m+k];
        }
      D[i*m+i] = 1 - sum;
      }
    
    // solve residual twist
    if( p > m-1 )
      {
      Xci = Xc;
      for(int i = 0; i < m; i++, Xci += p )
        for(int k = 0; k < p; k++ )
          if(sum_We[i] > 0 )
            Xci[k] /= sum_We[i];
      }

    // add vertex nudge
    Xci = Xc;
    for(int i = 0; i < m; i++, Xci += p )
      {
      double *Xk = X;
      for(int k = 0; k < m; k++, Xk += p )
        for(int h = 0; h < p; h++ )
          Xci[h] += D[i*m+k] * Xk[h];
      }

    if( prior_simplex )
      {
      Xci = Xc;
      double *X0i = X0;
      for(int i = 0; i < m; i++, Xci += p, X0i += p )
        {
        if( wX0[i] > 0 )
          {
          if( wX0[i] <= DBL_MAX ) // finite weight
            for(int h = 0; h < p; h++ )
              {
              Xci[h] = (Xci[h] + wX0[i]*X0i[h])/(1.0+wX0[i]);
              }
          else
            for(int h = 0; h < p; h++ )
              Xci[h] = X0i[h];
           }
        }
      }

    // copy Xc to X, and compute the norm of the difference
    double DX_F = 0, X_F = 0;
    for(int i = 0; i < m*p; i++ )
      {
      double f = X[i] - Xc[i];
      DX_F += f*f;
      X_F += X[i]*X[i];
      X[i] = Xc[i];
      }
    if( verbose > 0 )
      fprintf(stderr,"iter: %3d, |Delta X|_F/|X|_F = %g\n", r, sqrt(DX_F/X_F));

    if( DX_F/X_F < tolerance*tolerance ) break;
    } // for (int r = 0, ...)

  if(Beta)
    {
    double *beta_j = Beta, *Ycj = Yc;
    for(int j = 0; j < n; j++, beta_j += m, Ycj += p+1 )
      for(int i = 0; i < m; i++ )
        beta_j[i] = Ycj[i];
    }

  if(fitcone)  // reparameterize X and Beta
    {
    double *Xi = X + p;
    double *norm = Xc;
    for(int i = 1; i < m; i++, Xi += p )
      {
      norm[i] = 0;
      for(int k = 0; k < p; k++ )
        {
        Xi[k] -= X[k];  
        norm[i] += Xi[k]*Xi[k];
        }
      norm[i] = sqrt(norm[i]);
      for(int k = 0; k < p; k++ )
        Xi[k] /= norm[i];
      }
    double *beta_j = Beta;
    for(int j = 0; j < n; j++, beta_j += m )
      {
      beta_j[0] = 1.0;
      for(int i = 1; i < m; i++ )
        beta_j[i] *= norm[i];
      }
    }
  
  free(workspace);
  return 0;
}

// sfit() - older version of interface, without prior simplex (equivalent
// to sfit0() with prior simplex having all zero weights).
//

int
sfit (
  int n,                // data length
  int p,                // data dimensions
  const double *Y,      // data matrix, n x p
  const double *w,      // observation weight, n
  int m,                // simplex order
  int autoinit,         // 1: call initial_simplex, 0: use caller's X
  double *X,            // initial simplex (overwritten by the final), m x p
  double *Beta,         // affine combination coefficients, n x m

  double lambda,        // vertex assigment parameter
  double alpha,         // asymmetry weight
  int M_family,   // 0: least-squares, 1: Huber, 2: Biweight
  double gamma,         // robustness scale factor
                        // if < 0, use 1.345 for Huber and 4.685 for biweight
  int fitcone,          // if 1, treat first vertex as the apex (use separate
                        // robust scale estimate for the opposite facet).
                        // Otherwise, treat as simplex and pool scale est.
  
  int verbose,          // dump iteration steps to stderr
  double tolerance,     // convergence tolerance, ||\Delta X||F/||X||_F
  int iter_max,         // maximum iterations
  double Rtolerance     // smallest acceptable abs value of the diagonal entry
                        // of R still used for solving. 
  )
{
  double *X0 = NULL, *wX0 = NULL;

  return sfit0( n, p, Y, w, m, autoinit,
    X, wX0, X0, Beta, lambda, alpha, M_family, gamma,
    fitcone,verbose,tolerance, iter_max, Rtolerance);
}
// wrapper for indirection for calling from R
void
Rwrapper_sfit (
    int *n,
    int *p,
    double *Y,
    double *w,
    int *m,
    int *autoinit,
    double *lambda,
    double *alpha,
    int *M_family,
    double *gamma,
    int *fitcone,
    int *verbose,
    double *tolerance,
    int *iter_max,
    double *Rtolerance,
    double *X,          // output
    double *Beta        // output
    )
{
  sfit0(*n,*p, Y, w, *m, *autoinit, X, NULL, NULL, Beta,
      *lambda, *alpha, *M_family, *gamma, *fitcone,
      *verbose, *tolerance,
      *iter_max, *Rtolerance );
}

void
Rwrapper_sfit0 (
    int *n,
    int *p,
    double *Y,
    double *w,
    int *m,
    int *autoinit,
    double *lambda,
    double *alpha,
    int *M_family,
    double *gamma,
    int *fitcone,
    int *verbose,
    double *tolerance,
    int *iter_max,
    double *Rtolerance,
    double *X,          // output
    double *Beta,        // output
    double *wX0,
    double *X0
    )
{
  if (wX0 && !isfinite(wX0[0])) wX0 = NULL;
  if (X0 && isnan(X0[0])) X0 = NULL;

  sfit0(*n,*p, Y, w, *m, *autoinit, X, wX0, X0, Beta,
      *lambda, *alpha, *M_family, *gamma, *fitcone,
      *verbose, *tolerance,
      *iter_max, *Rtolerance );
}
#if CLI

#include <stdio.h>
#include <getopt.h>
#include "nio.h"

void
usage(char *argv0)
{
  fprintf(stderr,
  "usage:\n"
  " %s [options] [ <input> [ <output> ] ]\n"
  "\n"
  "-h, --help           : help\n"
  "\n",
  argv0 );
  exit(1);
}

int
main(int argc, char* argv[])
{
  /* option variable definitions */
  char *infile = "-";

  double lambda = 2;
  int M_family = 2;
  double alpha = 0.005;
  double gamma = -1; // default (huber 1.345, biweight 4.685)
  int fitcone = 1;   
  
  double convergence_tolerance = 1e-3;
  int maxit = 60;
  double Rtolerance = 1e-7;

  for(;;)
    {
    static struct option option_def[] = 
      {
        {"help", 0, 0, 'h'},
        { 0, 0, 0, 0 }
      };
    int c = getopt_long ( argc, argv, "h", option_def, 0 );
    if( c == -1 ) break;
    switch(c)
      {
      case 'h':
      default:
        usage(argv[0]);
        break;
      }
    }
  int m = atoi(argv[optind]);
  if(argc - optind > 1) infile = argv[optind+1];
   
  int n, p;
  double *y = get(infile, &n, &p);
  double *w = 0;
  double *b = (double*)malloc(sizeof(double)*m*n);
  double *x = (double*)malloc(sizeof(double)*m*p);
  
  int autoinit = 1;
  sfit(n, p, y, w, m, autoinit, x, b, lambda, alpha, M_family, gamma, fitcone,
      1, convergence_tolerance, maxit, Rtolerance );
  put("-", m, p, x );

  exit(0);
}
#endif
