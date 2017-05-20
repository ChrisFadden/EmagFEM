/*!\file

 *	Functions for Spectral Element Methods
 */

#ifndef SPEC1D_H
#define SPEC1D_H
#include "Constants.h"
#include "matrixIO.h"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

/*!
 * \brief Calculate points at the zeros of the polynomial
 * \param INPUT: N The order of the polynomial
 * \param OUTPUT:	xi	The interpolation points
 */
void calcLobattoPoints(const std::size_t N, Vec &xi) {
  // wi not actually used

  xi.resize(N);

  switch (N) {
  case 1:
    xi(0) = 0.0;
    // wi(0) = 4 / 3;
    return;
  case 2:
    xi(0) = -1.0 / std::sqrt(5.0);
    xi(1) = -xi(0);
    // wi(0) = 5/6;
    // wi(1) = wi(0);
    return;
  case 3:
    xi(0) = -std::sqrt(3.0 / 7.0);
    xi(1) = 0;
    xi(2) = -xi(0);
    // wi(0) = 49/90;
    // wi(1) = 32/45;
    // wi(2) = w(0);
    return;
  case 4:
    xi(0) = -0.76505532392946;
    xi(1) = -0.28523151648064;
    xi(2) = -xi(1);
    xi(3) = -xi(0);
    // wi(0) = 0.37847495629785;
    // wi(1) = 0.55485837703549;
    // wi(2) = wi(1);
    // wi(3) = wi(0);
    return;
  case 5:
    xi(0) = -0.83022389627857;
    xi(1) = -0.46884879347071;
    xi(2) = 0;
    xi(3) = -xi(1);
    xi(4) = -xi(0);
    // wi(0) = 0.27682604736157;
    // wi(1) = 0.43174538120986;
    // wi(2) = 0.48761904761905;
    // wi(3) = wi(1);
    // wi(4) = wi(0);
    return;
  case 6:
    xi(0) = -0.87174014850961;
    xi(1) = -0.59170018143314;
    xi(2) = -0.20929921790248;
    xi(3) = -xi(2);
    xi(4) = -xi(1);
    xi(5) = -xi(0);
    // wi(0) = 0.21070422714350;
    // wi(1) = 0.34112269248350;
    // wi(2) = 0.41245879465870;
    // wi(3) = wi(2);
    // wi(4) = wi(1);
    // wi(5) = wi(0);
    return;
  }
}

/*!
 * \brief 1-D linear Meshing
 * \param INPUT: x1 starting node
 * \param INPUT: x2 ending node
 * \param INPUT: N	Length of interval
 * \param INPUT: ratio of the last element / first element
 * \param OUTPUT: xe	array returning the node vector
 */
void elm_line1(const real x1, const real x2, const std::size_t N,
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
 * \brief Discretize the given domain using spectral elements
 * \param INPUT:	x0,x1 endpoints of the interval
 * \param INPUT:	Ne		Number of elements in the interval
 * \param INPUT:	ratio	Ratio of distances between nodes
 * \param INPUT:	np		Order of polynomial for each element
 * \param OUTPUT: xe		Element end nodes
 * \param OUTPUT:	xen		Element x interpolation nodes
 * \param OUTPUT: xien  Element xi interpolation nodes
 * \param OUTPUT:	xg		Global nodes
 * \param OUTPUT:	Cg		GLobal connectivity matrix
 * \param OUTPUT:	ng		number of unique global nodes
 */
void discretizeDomain(const real x0, const real x1, const std::size_t Ne,
		      const real ratio, const Vec &np, Vec &xe, Mat &xen,
		      Mat &xien, Vec &xg, Mat &Cg, std::size_t &ng) {

  elm_line1(x0, x1, Ne, ratio, xe);
  std::size_t m; // poly order
  Vec ti;	// tmp
  Vec xix;       // tmp
  Vec xi;

  xen.resize(Ne, 8);  // HARDCODED
  xien.resize(Ne, 8); // CHANGE THIS
  xen = 0 * xen;
  xien = 0 * xien;

  for (loop ll = 0; ll != Ne; ++ll) {
    m = np(ll);
    xi.resize(m + 1);
    xix.resize(m + 1);
    if (m > 1) {
      calcLobattoPoints(m - 1, ti);
      for (loop jj = 1; jj != m; ++jj) {
	xi(jj) = ti(jj - 1);
      }
    }
    xi(0) = -1.0;
    xi(m) = 1.0;

    for (loop jj = 0; jj != m + 1; ++jj) {
      xix(jj) =
	  0.5 * (xe(ll + 1) + xe(ll)) + 0.5 * (xe(ll + 1) - xe(ll)) * xi(jj);
    }

    for (loop jj = 0; jj != m + 1; ++jj) {
      xien(ll, jj) = xi(jj);
      xen(ll, jj) = xix(jj);
    }
  }

  ng = 1;
  Cg.resize(Ne, 8);
  Cg = 0 * Cg;
  xg.resize(8 * Ne);

  for (loop ll = 0; ll != Ne; ++ll) {
    ng = ng - 1;
    for (loop jj = 0; jj != np(ll) + 1; ++jj) {
      Cg(ll, jj) = ng;
      xg(ng) = xen(ll, jj);
      ng++;
    }
  }

  ng--;

  return;
}
/*!
 * \brief Create Vandermonde matrix for lobatto expansion
 * \param INPUT:  m		Order of the expansion
 * \param INPUT:	xi	node values for the given order of the expansion
 * \param OUTPUT:	Inverted Vandermonde matrix
 */
void vdm_modal_lob(const std::size_t m, const Vec xi, Mat &vdm) {
  vdm.resize(m + 1, m + 1);
  vdm = 0 * vdm;
  vdm(0, 0) = 1;
  vdm(m, m) = 1;

  real x;
  Mat vdmInv;

  for (loop jj = 1; jj != m + 1; ++jj) {
    x = xi(jj);
    vdm(0, jj) = (1 - x) / 2;
    vdm(1, jj) = 1 - x * x;
    vdm(m, jj) = (1 + x) / 2;

    if (m >= 3)
      vdm(2, jj) = 3 * (1 - x * x) * x;
    if (m >= 4)
      vdm(3, jj) = 3 * (1 - x * x) * (5 * x * x - 1) / 2;
    if (m >= 5)
      vdm(4, jj) = 5 * (1 - x * x) * (7 * x * x - 3) * x / 2;
    if (m >= 6)
      vdm(5, jj) =
	  15 * (1 - x * x) * (21 * std::pow(x, 4) - 14 * x * x + 1) / 8;
    if (m >= 7)
      vdm(6, jj) =
	  (1 - x * x) * (693 * std::pow(x, 4) - 630 * x * x + 105) * x / 8;
    if (m >= 8)
      vdm(7, jj) = (1 - x * x) * (3003 * std::pow(x, 6) -
				  3465 * std::pow(x, 4) + 945 * x * x - 35) /
		   16;
  }
  vdm = vdm.inverse();
}

/*!
 * \brief Create Element Diffusion matrix for lobatto expansion
 * \param INPUT:  m		Order of the expansion
 * \param INPUT:	dx  Element spacing
 * \param OUTPUT:	Element Diffusion matrix
 */
void edm_modal_lob(const std::size_t m, const real dx, Mat &edm) {
  edm.resize(m + 1, m + 1);
  edm *= 0;

  for (loop ii = 1; ii != m; ++ii) {
    edm(ii, ii) = 4 * ii * ii * (ii + 1) * (ii + 1) / ((2 * (ii + 1) - 1) * dx);
  }
  edm(0, 0) = 1 / dx;
  edm(0, m) = -1 / dx;
  edm(m, 0) = -1 / dx;
  edm(m, m) = 1 / dx;
}

/*!
 * \brief Create Element Mass matrix for lobatto expansion
 * \param INPUT:  m		Order of the expansion
 * \param INPUT:	dx  Element spacing
 * \param OUTPUT:	Element Mass matrix
 */
void emm_modal_lob(const std::size_t m, const real dx, Mat &emm) {
  emm.resize(m + 1, m + 1);
  emm *= 0;

  for (loop ii = 1; ii != m; ++ii) {
    emm(ii, ii) =
	dx * 2 * ii * ii * (ii + 1) * (ii + 1) /
	((2 * (ii + 1) - 3) * (2 * (ii + 1) - 1) * (2 * (ii + 1) + 1));
  }

  for (loop ii = 1; ii != m - 1; ++ii) {
    emm(ii, ii + 2) =
	-dx * ii * (ii + 1) * (ii + 2) * (ii + 3) /
	((2 * (ii + 1) - 1) * (2 * (ii + 1) + 1) * (2 * (ii + 1) + 3));
  }

  for (loop ii = 2; ii != m; ++ii) {
    emm(ii, ii - 2) =
	-dx * (ii - 2) * (ii - 1) * (ii) * (ii + 1) /
	((2 * (ii + 1) - 5) * (2 * (ii + 1) - 3) * (2 * (ii + 1) - 1));
  }

  emm(0, 0) = dx / 3;
  emm(0, 1) = dx / 3;
  emm(1, 0) = dx / 3;
  emm(0, 2) = -dx / 5;
  emm(2, 0) = -dx / 5;

  emm(0, m) = dx / 6;
  emm(1, m) = dx / 3;
  emm(2, m) = dx / 5;
  emm(m, 0) = dx / 6;
  emm(m, 1) = dx / 3;
  emm(m, 2) = dx / 5;
  emm(m, m) = dx / 3;
}

// INPUT should actually be after the parameter i.e.
// \param ne  INPUT:	blah blah blah

/*!
 * \brief Discretize the given domain using spectral elements
 * \param INPUT:	ne		Number of elements in the interval
 * \param INPUT: xe		Element end nodes
 * \param INPUT:	np		Order of polynomial for each element
 * \param INPUT:	ng		number of unique global nodes
 * \param INPUT:	Cg		Global connectivity matrix
 * \param INPUT: xien  Element xi interpolation nodes
 * \param INPUT: q0		 Boundary Condition
 * \param INPUT: fL		 Boundary Condition
 * \param INPUT: k		 Boundary Condition
 * \param INPUT: Src	 Source Excitation
 * \param OUTPUT:	A		 Matrix to be inverted
 * \param OUTPUT: b		 RHS of the Linear System
 */
void assembleSystem(const std::size_t ne, const Vec xe, const Vec np,
		    std::size_t ng, const Mat Cg, const Mat xien, const real q0,
		    const real fL, const real k, const Vec src, Mat &A,
		    Vec &b) {
  A.resize(ng + 1, ng + 1);
  b.resize(ng + 1);
  A *= 0;
  b *= 0;

  b(0) = q0 / k;

  loop m;
  Mat VdmInv;
  Mat Edm;
  Mat Emm;
  Vec xi;
  for (loop ll = 0; ll < ne; ++ll) {
    m = np(ll);
    xi.resize(m + 1);
    for (loop jj = 0; jj != m + 1; ++jj) {
      xi(jj) = xien(ll, jj);
    }

    vdm_modal_lob(m, xi, VdmInv);
    edm_modal_lob(m, xe(ll + 1) - xe(ll), Edm);
    emm_modal_lob(m, xe(ll + 1) - xe(ll), Emm);

    Edm = VdmInv * Edm * VdmInv.transpose();
    Emm = VdmInv * Emm * VdmInv.transpose();

    // Assemble FEM matrices
    for (loop ii = 0; ii != m + 1; ++ii) {
      for (loop jj = 0; jj != m + 1; ++jj) {
	A(Cg(ll, ii), Cg(ll, jj)) += Edm(ii, jj);
	b(Cg(ll, ii)) += Emm(ii, jj) * src(Cg(ll, jj)) / k;
      }
    }
  }

  // Boundary node
  m = np(ne);
  loop idx;
  for (loop ii = 0; ii != m; ++ii) {
    idx = Cg(ne - 1, ii);
    b(idx) -= Edm(ii, m) * fL;
  }
}

void test_Spec1D() {

  real x0 = 0;
  real x1 = 1;
  real k = 1.0;
  real fL = 0;
  real q0 = -1.0;
  real ratio = 1;

  std::size_t Ne = 4; // number of elements
  Vec np;	     // order of polynomial
  np.resize(Ne + 1);
  np(0) = 6;
  np(1) = 6;
  np(2) = 4;
  np(3) = 4;
  np(4) = 4;

  // Discretize domain
  Vec xe;
  Mat xen;
  Mat xien;
  Vec xg;
  Mat Cg;
  std::size_t ng;
  discretizeDomain(x0, x1, Ne, ratio, np, xe, xen, xien, xg, Cg, ng);

  // Set source
  Vec src;
  calcSrc(ng, src, xg);

  // Assemble System
  Mat A;
  Vec b;
  assembleSystem(Ne, xe, np, ng, Cg, xien, q0, fL, k, src, A, b);

  A.conservativeResize(ng, ng);
  b.conservativeResize(ng);
  // A.transposeInPlace();

  Vec x = A.llt().solve(b);

  x.conservativeResize(x.size() + 1);
  x(ng) = fL;
  DEBUG(x)

  return;
}

#endif
