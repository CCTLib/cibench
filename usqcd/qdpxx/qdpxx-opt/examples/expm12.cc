// $Id: expm12.cc,v 1.2 2004-07-27 05:38:37 edwards Exp $
/*! \file
 *  \brief 12-th order exponentiation of a lattice color matrix
 */

#include "qdp.h"
#include "examples.h"

using namespace QDP;

//! 12-th order exponentiation of a lattice color matrix
/*!
 * \ingroup gauge
 *
 *  In place  a_ = 1 + a_ + (1/2)*a_^2 + ...+ (1/n!)*(a_)^n  n = power
 *
 * Arguments:
 *
 *  \param a        LatticeColorMatrix          (Modify)
 */

void expm12(LatticeColorMatrix& a)
{
  START_CODE();
        
  // aux1 = aux2 = a;  a = ONE + a 
  LatticeColorMatrix aux1 = a;
  LatticeColorMatrix aux2 = a;
  LatticeColorMatrix aux3;
  a += 1;
  
  // Do a 12th order exponentiation
  for(int i= 2; i <= 12; ++i)
  {
    Real dummy = Real(1)/Real(i);

    aux3 = aux2 * aux1;
    aux2 = aux3 * dummy;
    a += aux2;
  }
  
  END_CODE();
}
