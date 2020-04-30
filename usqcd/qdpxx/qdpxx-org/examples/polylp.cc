// $Id: polylp.cc,v 1.2 2004-07-27 05:38:37 edwards Exp $
/*! \file
 *  \brief Calculate the global normalized sum of the Polyakov loop
 */

#include "qdp.h"
#include "examples.h"

using namespace QDP;

//! Return the value of the average plaquette normalized to 1
/*!
 * \ingroup glue
 *
 * \param u          gauge field (Read)
 * \param poly_loop  Polyakov loop average in direction mu (Write) 
 * \param mu         direction of Polyakov loop (Read)
 */

void polylp(const multi1d<LatticeColorMatrix>& u, DComplex& poly_loop, int mu)
{
  START_CODE();
        
  // Initial Polyakov loop
  LatticeColorMatrix poly = u[mu];

  for(int n = 1; n < Layout::lattSize()[mu]; ++n)  // run over all links in mu dir
  {    
    LatticeColorMatrix tmp = shift(poly, FORWARD, mu);
    poly = u[mu] * tmp;
  }

  /* Take the trace and sum up */
  poly_loop = sum(trace(poly)) / Double(Nc*Layout::vol());

  END_CODE();
}
