//  $Id: mesons_w.cc,v 1.12 2007-02-22 15:58:30 bjoo Exp $
/*! \file
 *  \brief Meson 2-pt functions
 */

#include "examples.h"

using namespace QDP;

//! Function object used for constructing the time-slice set
class TimeSliceFunc : public SetFunc
{
public:
  TimeSliceFunc(int dir): dir_decay(dir) {}

  int operator() (const multi1d<int>& coordinate) const {return coordinate[dir_decay];}
  int numSubsets() const {return Layout::lattSize()[dir_decay];}

  int dir_decay;

private:
  TimeSliceFunc() {}  // hide default constructor
};

 
//! Meson 2-pt functions
/* This routine is specific to Wilson fermions!
 *
 * Construct meson propagators
 * The two propagators can be identical or different.
 *
 * \param quark_prop_1 -- first quark propagator ( Read )
 * \param quark_prop_2 -- second (anti-) quark propagator ( Read )
 * \param meson_propagator -- Ns^2 mesons ( Modify )
 * \param t_source -- cartesian coordinates of the source ( Read )
 * \param j_decay -- direction of the exponential decay ( Read )
 *
 *        ____
 *        \
 * m(t) =  >  < m(t_source, 0) m(t + t_source, x) >
 *        /
 *        ----
 *          x
 */

void mesons(const LatticePropagator& quark_prop_1, const LatticePropagator& quark_prop_2, 
	    multi1d< multi1d<Real> >& meson_propagator, 
	    const multi1d<int>& t_source, int j_decay)
{
  START_CODE();

  // Create the time-slice set
  Set timeslice;
  timeslice.make(TimeSliceFunc(j_decay));

  // Length of lattice in j_decay direction
  int length = timeslice.numSubsets();

  int t0 = t_source[j_decay];
  int G5 = Ns*Ns-1;
  multi1d<Double> hsum(length);

  meson_propagator.resize(Ns*Ns);  

  // Contruct the antiquark prop
  LatticePropagator anti_quark_prop =  Gamma(G5) * quark_prop_2 * Gamma(G5);

  for(int n = 0; n < (Ns*Ns); ++n)
  {
    // Initialize the propagator so that we just add to it below
    meson_propagator[n].resize(length);
    meson_propagator[n] = zero;

    LatticeReal psi_sq = real(trace(adj(anti_quark_prop) * Gamma(n) * quark_prop_1 * Gamma(n)));

    // Do a slice-wise sum.
    hsum = sumMulti(psi_sq, timeslice);

    for(int t = 0; t < length; ++t)
    {
      int t_eff = (t - t0 + length) % length;

      meson_propagator[n][t_eff] += Real(hsum[t]);
    }
  }

  END_CODE();
}
