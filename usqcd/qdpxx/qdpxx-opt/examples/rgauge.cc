// $Id: rgauge.cc,v 1.2 2004-07-27 05:38:37 edwards Exp $
/*! \file
 *  \brief Do a random gauge transformation on the u fields
 */

#include "qdp.h"
#include "examples.h"

using namespace QDP;

//! Do a random gauge transformation on the u fields
/*!
 * \ingroup gauge
 *
 * Arguments:
 *
 *  \param u          Gauge field                   (Modify)
 *  \param g          Gauge transformation matrices (Write)
 */

void rgauge(multi1d<LatticeColorMatrix>& u, LatticeColorMatrix& g)
{
  START_CODE();
  
  // g = exp(A), where A = random traceless antihermitian matrix.
  
  // NOTE!!!!: we are using a hack. The twelth_order exponentiation followed
  //   by a reunitarization should be replace by an exact exponentiation.
  //   So, the gauge field distriubution is NOT from Haar measure, but
  //   if someone can think of a really good reason why it should be Haar
  //   measure, then something can be done about it.
  
  gaussian(g);

  taproj(g);
  expm12(g);
  reunit(g);
    
  for(int mu = 0; mu < Nd; ++mu)
  {
    LatticeColorMatrix u_tmp = g * u[mu];
    u[mu] = u_tmp * shift(adj(g), FORWARD, mu);
  }
    
  END_CODE();
}
