/*!\file

 *	Functions for FEM in 1D from Finite and Spectral Element Methods with
 *MATLAB
 */

#ifndef FEM1D_H
#define FEM1D_H
#include "Constants.h"
#include "matrixIO.h"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
/*!
 * \brief 1-D linear Meshing
 * \param INPUT: x1 starting node
 * \param INPUT: x2 ending node
 * \param INPUT: N	Length of interval
 * \param INPUT: ratio of the last element / first element
 * \param OUTPUT: xe	array returning the node vector
 */
void elm_line1(const node x1, const node x2, const std::size_t N,
	       const real ratio, Vec &xe) {

  xe.resize(N + 1);

  // one element
  if (N == 1) {
    xe(0) = x1;
    xe(1) = x2;
    return;
  }

  real alpha, factor;
  if (ratio == 1) {
    alpha = 1.0;
    factor = 1.0 / N;
  } else {
    alpha = std::pow(ratio, 1.0 / (N - 1));
    factor = (1.0 - alpha) / (1.0 - std::pow(alpha, N));
  }

  real deltax = (x2 - x1) * factor;
  xe(0) = x1;
  for (loop ii = 1; ii < N + 1; ++ii) {
    xe(ii) = xe(ii - 1) + deltax;
    deltax = deltax * alpha;
  }

  return;
}

/*!
 * \brief Sets the source function values
 * \param INPUT:	N Length of Interval
 * \param OUTPUT:	src Vector
 */
void calcSrc(const std::size_t N, Vec &src, const Vec &xe) {
  src.resize(N + 1);

  for (loop ii = 0; ii < N + 1; ++ii) {
    src(ii) = 10.0 *
	      std::exp(-5.0 * std::pow(xe(ii), 2) / std::pow(xe(N) - xe(0), 2));
  }
}

/*!
 * \brief Assembles the Finite Element Matrix A and RHS b for the steady state
 *	diffusion (Poisson) equation \f$k\frac{d^2f}{dx^2} + s(x)=0\f$
 *        with Neumann BC
 *				\f$q_0 = -k \frac{df}{dx}[x=0]\f$ and Diriclet
 *BC
 *\f$f(x=L) = f_L\f$
 * \param INPUT:	N	Length of Interval
 * \param INPUT:	xe	node listing
 * \param INPUT:	q0	BC parameter
 * \param INPUT:	fL	BC parameter
 * \param INPUT:	k		BC parameter
 * \param INPUT:	s		Src Vector
 * \param OUTPUT:	A		FEM Matrix
 * \param OUTPUT:	b		FEM Src
 */
void assembleSystem_Diff(const std::size_t N, const Vec &xe, const real q0,
			 const real fL, const real k, const Vec &s, Mat &A,
			 Vec &b) {

  b.resize(N);
  b = 0 * b;
  b(0) = q0 / k; // Neumann BC

  A.resize(N, N);
  A = 0 * A;

  // Handle start node
  real A11 = 1 / (xe(1) - xe(0));
  real A12 = -A11;
  real A21 = A12;
  real A22 = A11;
  real B11 = (xe(1) - xe(0)) / 3.0;
  real B12 = 0.5 * B11;
  real B21 = B12;
  real B22 = B11;

  A(0, 0) += A11;
  A(0, 1) += A12;
  A(1, 1) += A22;

  b(0) += (B11 * s(0) + B12 * s(1)) / k;
  b(1) += (B21 * s(0) + B22 * s(1)) / k;

  // Loop over nodes
  for (loop ii = 1; ii < N - 1; ++ii) {
    A11 = 1 / (xe(ii + 1) - xe(ii));
    A12 = -A11;
    A21 = A12;
    A22 = A11;

    B11 = (xe(ii + 1) - xe(ii)) / 3.0;
    B12 = 0.5 * B11;
    B21 = B12;
    B22 = B11;

    A(ii, ii) += A11;
    A(ii, ii + 1) += A12;
    A(ii, ii - 1) += A21;
    A(ii + 1, ii + 1) += A22;

    b(ii) += (B11 * s(ii) + B12 * s(ii + 1)) / k;
    b(ii + 1) += (B21 * s(ii) + B22 * s(ii + 1)) / k;
  }

  //	 Handle End Node
  A11 = 1 / (xe(N) - xe(N - 1));
  B11 = (xe(N) - xe(N - 1)) / 3.0;
  B12 = 0.5 * B11;
  A12 = -A11;
  A21 = A12;
  A(N - 1, N - 2) += A21;
  A(N - 1, N - 1) += A11;
  b(N - 1) += (B11 * s(N - 1) + B12 * s(N)) / k - A12 * fL;
}

/*!
 * \brief Finite and Spectral Elements w/ MATLAB pg. 31 Poisson Eqn
 */
void test_fem1d_Diff() {

  Vec xe;  // Grid Nodes
  Vec src; // Source Vector
  Mat A;
  Vec b;

  // Simulation Parameters
  node x0 = 0;
  node x1 = 1;
  std::size_t N = 8;
  real ratio = 1;
  real k = 1.0;
  real q0 = -1.0;
  real fL = 0;

  // Grid elements
  elm_line1(x0, x1, N, ratio, xe);

  // Source values
  calcSrc(N, src, xe);

  // Setup Linear System
  assembleSystem_Diff(N, xe, q0, fL, k, src, A, b);

  // Solve System
  Vec x = A.colPivHouseholderQr().solve(b);
  // vec2bin("test.dat", x);
  return;
}

/*!
 * \brief Assembles the Finite Element Matrix A and RHS b for the
 *	Helmholtz equation \f$k\frac{d^2f}{dx^2} + \alpha f =0\f$
 *        with Neumann BC
 *				\f$q_0 = -k \frac{df}{dx}[x=0]\f$ and Diriclet
 *BC
 *\f$f(x=L) = f_L\f$
 * \param INPUT:	N	Length of Interval
 * \param INPUT:	xe	node listing
 * \param INPUT:	q0	BC parameter
 * \param INPUT:	fL	BC parameter
 * \param INPUT:	k		BC parameter
 * \param INPUT: \f$\alpha\f$ wavenumber squared
 * \param OUTPUT:	A		FEM Matrix
 * \param OUTPUT:	b		FEM Src
 */
void assembleSystem_Helm(const std::size_t N, const Vec &xe, const real q0,
			 const real fL, const real k, const real alpha, Mat &A,
			 Vec &b) {

  b.resize(N);
  b = 0 * b;
  b(0) = q0 / k; // Neumann BC

  A.resize(N, N);
  A = 0 * A;

  // Handle start node
  real A11 = 1 / (xe(1) - xe(0));
  real A12 = -A11;
  real A21 = A12;
  real A22 = A11;
  real B11 = (xe(1) - xe(0)) / 3.0;
  real B12 = 0.5 * B11;
  real B21 = B12;
  real B22 = B11;

  A(0, 0) += A11 - alpha * B11;
  A(0, 1) += A12 - alpha * B12;
  A(1, 1) += A22 - alpha * B22;

  // Loop over nodes
  for (loop ii = 1; ii < N - 1; ++ii) {
    A11 = 1 / (xe(ii + 1) - xe(ii));
    A12 = -A11;
    A21 = A12;
    A22 = A11;

    B11 = (xe(ii + 1) - xe(ii)) / 3.0;
    B12 = 0.5 * B11;
    B21 = B12;
    B22 = B11;

    A(ii, ii) += A11 - alpha * B11;
    A(ii, ii + 1) += A12 - alpha * B12;
    A(ii, ii - 1) += A21 - alpha * B21;
    A(ii + 1, ii + 1) += A22 - alpha * B22;
  }

  //	 Handle End Node
  A11 = 1 / (xe(N) - xe(N - 1));
  B11 = (xe(N) - xe(N - 1)) / 3.0;
  B12 = 0.5 * B11;
  A12 = -A11;
  A21 = A12;

  A(N - 1, N - 2) += A21 - alpha * B21;
  A(N - 1, N - 1) += A11 - alpha * B11;
  b(N - 1) += b(N - 1) - (A12 - alpha * B12) * fL;
}

/*!
 * \brief Finite and Spectral Elements w/ MATLAB pg. 50 Helmholtz Eqn
 */
void test_fem1d_Helm() {
  Vec xe;  // Grid Nodes
  Vec src; // Source Vector
  Mat A;
  Vec b;

  // Simulation Parameters
  node x0 = 0;
  node x1 = 1;
  std::size_t N = 20;
  real ratio = 1.0;
  real k = 1.0;
  real q0 = -1.0;
  real fL = -0.2;
  real alpha = 87.4;

  elm_line1(x0, x1, N, ratio, xe);

  // Setup Linear System
  assembleSystem_Helm(N, xe, q0, fL, k, alpha, A, b);

  // Solve System
  Vec x = A.colPivHouseholderQr().solve(b);
  std::cout << x << std::endl;
}

/*!
 * \brief Assembles the Finite Element Matrix A and RHS b with Quadratic
 *Elements
 * for the steady state
 *	diffusion (Poisson) equation \f$k\frac{d^2f}{dx^2} + s(x)=0\f$
 *        with Neumann BC
 *				\f$q_0 = -k \frac{df}{dx}[x=0]\f$ and Diriclet
 *BC
 *\f$f(x=L) = f_L\f$
 * \param INPUT:	N	Length of Interval
 * \param INPUT:	xe	node listing
 * \param INPUT:	\f$\Beta\f$	Location of center node (0 = midpoint)
 * \param INPUT:	q0	BC parameter
 * \param INPUT:	fL	BC parameter
 * \param INPUT:	k		BC parameter
 * \param INPUT:	s		Src Vector
 * \param OUTPUT:	A		FEM Matrix
 * \param OUTPUT:	b		FEM Src
 */
void assembleSystem_QuadElm_Diff(const std::size_t N, const Vec &xe,
				 const real beta, const real q0, const real fL,
				 const real k, const Vec &s, Mat &A, Vec &b) {

  b.resize(2 * N + 1);
  b = 0 * b;
  b(0) = q0 / k; // Neumann BC

  A.resize(2 * N + 1, 2 * N + 1);
  A = 0 * A;

  Eigen::Matrix<real, 3, 3> Aii;
  Eigen::Matrix<real, 3, 3> AiiConst;
  Eigen::Matrix<real, 3, 3> Bii;
  Eigen::Matrix<real, 3, 3> BiiConst;
  Eigen::Matrix<real, 4, 1> xi;
  Eigen::Matrix<real, 4, 1> wi;
  Eigen::Matrix<real, 3, 1> psi;

  real cf;
  loop ii1, ii2, ii3;

  // Precompute matrix elements
  real beta2pos = std::pow(1 + beta, 2);
  real beta2neg = std::pow(1 - beta, 2);
  real beta2 = std::pow(beta, 2);
  AiiConst(0, 0) = (4 + 3 * beta2pos) / beta2pos;
  AiiConst(0, 1) = -8 / ((1 + beta) * (1 - beta2));
  AiiConst(0, 2) = (4 - 3 * (1 - beta2)) / (1 - beta2);
  AiiConst(1, 0) = AiiConst(0, 1);
  AiiConst(1, 1) = 16 / std::pow((1 - beta2), 2);
  AiiConst(1, 2) = -8 / ((1 - beta) * (1 - beta2));
  AiiConst(2, 0) = AiiConst(0, 2);
  AiiConst(2, 1) = AiiConst(1, 2);
  AiiConst(2, 2) = (4 + 3 * beta2neg) / beta2neg;

  // nodes
  xi(0) = -1.0;
  xi(1) = -1 / std::sqrt(5);
  xi(2) = -xi(1);
  xi(3) = -xi(0);

  // weights
  wi(0) = 1.0 / 6.0;
  wi(1) = 5.0 / 6.0;
  wi(2) = wi(1);
  wi(3) = wi(0);

  // initialize Bii
  for (loop jj = 0; jj < 3; ++jj) {
    for (loop kk = 0; kk < 3; ++kk) {
      BiiConst(jj, kk) = 0;
    }
  }

  for (loop qq = 0; qq < 4; ++qq) {
    psi(0) = (xi(qq) - beta) * (xi(qq) - 1) / (2 * (1 + beta));
    psi(1) = (1 - std::pow(xi(qq), 2)) / (1 - std::pow(beta, 2));
    psi(2) = (xi(qq) - beta) * (xi(qq) + 1) / (2 * (1 - beta));

    BiiConst(0, 0) += psi(0) * psi(0) * wi(qq);
    BiiConst(0, 1) += psi(0) * psi(1) * wi(qq);
    BiiConst(0, 2) += psi(0) * psi(2) * wi(qq);
    BiiConst(1, 1) += psi(1) * psi(1) * wi(qq);
    BiiConst(1, 2) += psi(1) * psi(2) * wi(qq);
    BiiConst(2, 2) += psi(2) * psi(2) * wi(qq);
  }

  BiiConst(1, 0) = BiiConst(0, 1);
  BiiConst(2, 0) = BiiConst(0, 2);
  BiiConst(2, 1) = BiiConst(1, 2);

  // Loop over nodes
  for (loop ii = 0; ii < N; ++ii) {
    cf = 1.0 / (3.0 * (xe(ii + 1) - xe(ii)));

    Bii = ((xe(ii + 1) - xe(ii)) / 2) * BiiConst;

    Aii = AiiConst * cf;

    // a = A(i,i)
    // b = A(i,i+1)
    // c = A(i,i+2)
    // d = A(i,i-1)
    // e = A(i,i-2)

    ii1 = 2 * ii;
    ii2 = 2 * ii + 1;
    ii3 = 2 * ii + 2;

    A(ii1, ii1) += Aii(0, 0);
    A(ii1, ii1 + 1) += Aii(0, 1);
    A(ii1, ii1 + 2) += Aii(0, 2);

    A(ii2, ii2 - 1) += Aii(1, 0);
    A(ii2, ii2) += Aii(1, 1);
    A(ii2, ii2 + 1) += Aii(1, 2);

    A(ii3, ii3 - 2) += Aii(2, 0);
    A(ii3, ii3 - 1) += Aii(2, 1);
    A(ii3, ii3) += Aii(2, 2);

    b(ii1) +=
	(Bii(0, 0) * s(ii1) + Bii(0, 1) * s(ii2) + Bii(0, 2) * s(ii3)) / k;
    b(ii2) +=
	(Bii(1, 0) * s(ii1) + Bii(1, 1) * s(ii2) + Bii(1, 2) * s(ii3)) / k;
    b(ii3) +=
	(Bii(2, 0) * s(ii1) + Bii(2, 1) * s(ii2) + Bii(2, 2) * s(ii3)) / k;

  } // loop over nodes

  // Handle boundary nodes
  b(2 * N - 2) -= A(2 * N - 2, 2 * N) * fL;
  b(2 * N - 1) -= A(2 * N - 1, 2 * N) * fL;
} // end Quad Elm function

/*!
 * \brief Finite and Spectral Elements w/ MATLAB pg. 79 Poisson Eqn w/ Quadratic
 * elements
 */

void test_fem1d_QuadElm_Diff() {

  Vec xe;
  Vec xg;
  Vec src;
  Mat A;
  Vec b;

  // Simulation Parameters
  node x0 = 0;
  node x1 = 1;
  real k = 1;
  real q0 = -1;
  real fL = 0;
  std::size_t N = 5;
  assert(N >= 5);
  std::size_t Ng = 2 * N;
  real ratio = 1.0;
  real beta = 0;

  elm_line1(x0, x1, N, ratio, xe);

  xg.resize(Ng + 1);
  loop jj = 0;
  for (loop ii = 0; ii < N; ii++) {
    xg(jj) = xe(ii);
    ++jj;
    xg(jj) = 0.5 * (xe(ii) + xe(ii + 1)) + 0.5 * (xe(ii + 1) - xe(ii)) * beta;
    ++jj;
  }

  xg(Ng) = xe(N);

  // Source values
  calcSrc(Ng, src, xg);

  assembleSystem_QuadElm_Diff(N, xe, beta, q0, fL, k, src, A, b);

  // Solve System
  Vec x = A.llt().solve(b);
  x(Ng) = fL;
  DEBUG(x)
  return;
}

#endif
