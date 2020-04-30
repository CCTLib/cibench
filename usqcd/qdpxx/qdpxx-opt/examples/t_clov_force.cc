// $Id: t_skeleton.cc,v 1.5 2003-06-07 19:09:32 edwards Exp $
/*! \file
 *  \brief Skeleton of a QDP main program
 */

#include "qdp.h"

#include "examples.h"

using namespace QDP;

void deriv_loops(const multi1d<LatticeColorMatrix>& u,
		 const int mu, const int nu, const int cb,
		 LatticeColorMatrix& ds_u_mu,
		 LatticeColorMatrix& ds_u_nu,
		 const LatticeColorMatrix& Lambda)
{
  
  // New thingie - now assume Lambda lives only on sites with checkerboard 
  // CB
  //            Lambda
  //   0           X           0           x = cb, O = 1-cb
  //
  //
  // Lambda                 Lambda 
  //   X           0           X
  //
  //
  //            Lambda 
  //   0           X           0
  //
  // So I can only construct 4 out of the 8 staples on the sites
  // that have CB and the OTHER 4 of the 8 staples on sites with 
  // 1-cb
  //
  
  // Sites with CB first:
  //
  
  LatticeColorMatrix staple_for;
  LatticeColorMatrix staple_back;
  LatticeColorMatrix staple_left;
  LatticeColorMatrix staple_right;
  
  LatticeColorMatrix u_nu_for_mu = shift(u[nu],FORWARD, mu); // Can reuse these later
  LatticeColorMatrix u_mu_for_nu = shift(u[mu],FORWARD, nu);
  LatticeColorMatrix Lambda_xplus_mu = shift(Lambda, FORWARD, mu);
  LatticeColorMatrix Lambda_xplus_nu = shift(Lambda, FORWARD, nu);
  LatticeColorMatrix Lambda_xplus_muplusnu = shift(Lambda_xplus_mu, FORWARD, nu);
  
  LatticeColorMatrix u_tmp3;
  
  LatticeColorMatrix ds_tmp_mu;
  LatticeColorMatrix ds_tmp_nu;
  {
    LatticeColorMatrix up_left_corner;
    LatticeColorMatrix up_right_corner;
    LatticeColorMatrix low_right_corner;
    LatticeColorMatrix low_left_corner;
    
    //   u_tmp1 =   <-------
    //              |
    //              |                       ON ALL CHECKERBOARDS
    //              |                       (Because it's used in staples)  
    //              V
    up_left_corner = adj(u_mu_for_nu)*adj(u[nu]);
    
    
    //
    //              <------^
    //                     |
    //                     |
    //                     |
    
    up_right_corner = u_nu_for_mu*adj(u_mu_for_nu);
    
    //                       |
    //                       |
    //                       |
    //                       V
    //                <------
    low_right_corner = adj(u_nu_for_mu)*adj(u[mu]);
    
    //
    //                    ^
    //  low left corner=  |                         ON ALL CHECKBERBOARDS
    //                    |                         (Because it's used in the staples)
    //                    |  
    //                     <-------
    low_left_corner = adj(u[mu])*u[nu];
    
    
    // Now compute the terms of the force:
    // 
    // Altogether 8 terms. 4 Upwards with + sign, and 4 Downwards with - sign
    //                     4 terms use staples and 4 don't
    
    // NON STAPLE TERMS FIRST:
    
    // 1) mu links
    //
    //    <-------  X (CB)                      <--------
    //    |         ^                           |
    //    |         |        re  use  u_tmp1 =  | 
    //    V         |                           V
    //   CB       1-CB       
    u_tmp3[rb[cb]] = u_nu_for_mu*Lambda_xplus_muplusnu;
    ds_u_mu[rb[cb]] = u_tmp3*up_left_corner;
    
    //    nu links
    //    X
    //     <------
    //     | 
    //     |
    //     |
    //     V-----> CB
    //   (1-CB)   
    //
    u_tmp3[rb[1-cb]] = adj(u_mu_for_nu)*Lambda_xplus_nu;
    
    // accumulate into ds_tmp_nu and shift everything together at the end
    ds_tmp_nu[rb[1-cb]] = u_tmp3*adj(low_left_corner);
    
    
    
    // 2)  mu links
    //
    //  CB
    //    X <------ 
    //    |        ^       re use u[nu](x+mu) = u_nu_for_mu
    //    |        |       re use u[mu](x+nu) = u_mu_for_nu
    //    V        |
    //    1-CB    CB        
    u_tmp3[rb[1-cb]] = Lambda_xplus_nu*adj(u[nu]);
    ds_u_mu[rb[1-cb]] = up_right_corner * u_tmp3;
    
    //      nu_links
    //    
    //     <------
    //     | 
    //     |
    //     |
    //   X V----->1-CB
    //   (CB)   
    //
    u_tmp3[rb[cb]] = up_left_corner*Lambda;
    //
    // accumulate into ds_tmp_nu and shift everything together at the end
    ds_tmp_nu[rb[cb]] = u_tmp3*u[mu];
    
    
    
    
    //
    // Terms 3) and 4)
    //
    // These last two can be done on the other checkerboard and then shifted together. at the very end...
    //
    //  CB      1-CB          
    //    ^       |   
    //    |       |     
    //    |       V
    //    <-------X CB
    
    
    // 3) Mu links
    //
    //  Compunte with low_left_corner:     ^           |
    //                                     |           | 
    //               low_left_corner    =  |           | 
    //                                     |           V
    //                              (1-CB) <--------   X CB
    u_tmp3[rb[1-cb]] = adj(u_nu_for_mu)*Lambda_xplus_mu;
    //
    // accumulate into ds_tmp_mu and shift at the end.
    ds_tmp_mu[rb[1-cb]] = u_tmp3*low_left_corner;
    
    // Nu links
    // 
    //  CB    ------>                                     ------->
    //                |                                          |
    //                |      reuse adj(up_right_corner):         |
    //                |                                          |
    //                V                                          V
    //  1-CB   <------X
    u_tmp3[rb[1-cb]] = adj(up_right_corner)*Lambda_xplus_mu;
    ds_u_nu[rb[1-cb]] = u_tmp3*adj(u[mu]);
    
    
    
    // 4) Mu links
    //
    //  1-CB      CB
    //   ^        |
    //   |        |        reuse = u[nu](x+mu) = u_nu_for_mu
    //   |        |
    //   X <----- V 1-CB
    //   CB
    u_tmp3[rb[cb]] = low_right_corner*Lambda;
    //
    // accumulate into ds_tmp_mu and shift at the end.
    ds_tmp_mu[rb[cb]] = u_tmp3*u[nu];
    
    
    
    // Nu links
    // 
    //  1-CB   ------> X                               
    //               |                                          |
    //               |       reuse low_right_corner:            |
    //               |                                          |
    //               V                                          V
    //   CB   <------                                    <------
    u_tmp3[rb[cb]] =    u_mu_for_nu*Lambda_xplus_muplusnu;
    ds_u_nu[rb[cb]] =   u_tmp3*low_right_corner;
    
    
    //  ds_tmp_mu now holds the last 2 terms, one on each of its checkerboards, but Now I need
    //  to shift them both together onto ds_u_mu
    //  I'll keep them in ds_tmp_mu right, bearing in mind I'll need to bring
    //  them in with a -ve contribution...
    
    
    // STAPLE TERMS:   
    
    // Construct the staples
    
    //  Staple_for =  <--------
    //                |       ^
    //                |       |             ON ALL CHECKERBOARDS
    //                |       |             
    //                V       |
    staple_for = u_nu_for_mu*up_left_corner;
    
    
    // Staple_right =   <-----             ON ALL CHECKERBOARDS
    //                 |
    //                 |
    //                 V
    //                 ----->
    staple_right = up_left_corner*u[mu];
    
    
    
    //                 ----->
    //                       |
    //                       |
    //                       |
    //                <----- V
    staple_left  = u_mu_for_nu*low_right_corner;
    
    
    
    
    //  Staple_back =  ^       |
    //                 |       |            ON ALL CHECKERBOARDS
    //                 |       |
    //                 <------ V
    //                      
    staple_back = adj(u_nu_for_mu)*low_left_corner;
    
  }  // Corner pieces go away here
  
  // 5) Mu links
  //
  //    <------- 
  //    |        ^
  //    |        |     use computed staple
  //    V        |
  //    x        
  //   CB       1-CB  
  ds_u_mu[rb[cb]] += staple_for*Lambda;
  
  
  //  Nu links
  //
  //     CB   <---- 1-CB
  //        |
  //        |                use staple_right
  //        V
  //    1-CB  -----> X CB
  //
  //  Accumulate into ds_tmp_nu and shift at the end.
  
  ds_tmp_nu[rb[1-cb]] += staple_right*Lambda_xplus_mu;
  
  
  // 6)  Mu links
  //
  //    <------- 
  //    |        ^
  //    |        |    re use computed staple 
  //    V        |
  //   1-CB      X CB	  
  
  ds_u_mu[rb[1-cb]] += Lambda_xplus_mu*staple_for;
  
  
  
  //  Nu links
  // 
  //      <----  X CB
  //     |
  //     |                     use adj(staple_right)
  //     |
  // CB  V ----> (1-CB)
  ds_tmp_nu[rb[cb]] += Lambda_xplus_muplusnu * staple_right;
  
  
  // 7) Mu links
  //
  //   CB      1-CB
  //  X         
  //    ^       |
  //    |       |   re use computed staple 
  //    |       |
  //    <-------V 
  //
  //  Accumulate into ds_tmp_mu and shift at the end.
  ds_tmp_mu[rb[1-cb]] += staple_back*Lambda_xplus_nu;
  
  //   Now for nu
  //
  //   (1-CB)  -----> CB         use adj(staple_left)
  //                |
  //                |
  //                V
  //      CB X <----
  //
  ds_u_nu[rb[cb]] += staple_left*Lambda;
  
  // 8) Mu links
  //
  //  1-CB      X CB
  //    ^       |
  //    |       |  reuse computed staple 
  //    |       |
  //    <-------V
  //
  // Accumulate into ds_tmp_mu and shift at the end
  ds_tmp_mu[rb[cb]] += Lambda_xplus_muplusnu * staple_back;
  
  // Now for Nu
  // 
  //    CB X ------> (1-CB)
  //               |
  //               |
  //               |
  //               V
  // 1-CB  <------- CB
  ds_u_nu[rb[1-cb]] += Lambda_xplus_nu * staple_left;
  
  // Now shift the accumulated pieces to mu and nu
  // 
  // Hope that this is not too slow as an expression
  ds_u_mu -= shift(ds_tmp_mu, BACKWARD, nu);
  ds_u_nu -= shift(ds_tmp_nu, BACKWARD, mu);   
}

void deriv(const multi1d<LatticeColorMatrix>& u,
	   multi1d<LatticeColorMatrix>& ds_u, 
	   const LatticeFermion& chi, const LatticeFermion& psi, 
	   int cb)
  {
    // Do I still need to do this?
    if( ds_u.size() != Nd ) { 
      ds_u.resize(Nd);
    }
    
    ds_u = zero;
    
    // Now compute the insertions
    for(int mu=0; mu < Nd; mu++) {
      for(int nu = mu+1; nu < Nd; nu++) {
	
	// These will be appropriately overwritten - no need to zero them.
	// Contributions to mu links from mu-nu clover piece
	LatticeColorMatrix ds_tmp_mu; 
	
	// -ve contribs  to the nu_links from the mu-nu clover piece 
	// -ve because of the exchange of gamma_mu gamma_nu <-> gamma_nu gamma_mu
	LatticeColorMatrix ds_tmp_nu;
	
	// The weight for the terms
	// I am going to assume for this test that the clover coeff is 1.
	Real factor = (Real(-1)/Real(8));

	// Get gamma_mu gamma_nu psi -- no saving here, from storing shifts because
	// I now only do every mu, nu pair only once.

	int mu_nu_index = (1 << mu) + (1 << nu); // 2^{mu} 2^{nu}
	LatticeFermion ferm_tmp = Gamma(mu_nu_index)*psi;
	LatticeColorMatrix s_xy_dag = traceSpin(outerProduct(ferm_tmp,chi));
	s_xy_dag *= Real(factor);

	// Compute contributions
	deriv_loops(u, mu, nu, cb, ds_tmp_mu, ds_tmp_nu, s_xy_dag);

	// Accumulate them
	ds_u[mu] += ds_tmp_mu;
	ds_u[nu] -= ds_tmp_nu;


      }
    }


    // Clear out the deriv on any fixed links
    //    (*this).getFermBC().zero(ds_u);
    // For this outlining, ignore this 
    
  }



int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the lattice size
  // NOTE: in general, the user will need/want to set the
  // lattice size at run-time
  multi1d<int> nrow(Nd);
  nrow[0] = 8;
  nrow[1] = 8; 
  nrow[2] = 8;
  nrow[3] = 8;

  // Insert the lattice size into the Layout
  // There can be additional calls here like setting the number of
  // processors to use in an SMP implementation
  Layout::setLattSize(nrow);

  // Create the lattice layout
  // Afterwards, QDP is useable
  Layout::create();

  // Do some wonderful and amazing things - impress your friends
  multi1d<LatticeColorMatrix> u(Nd);
  multi1d<LatticeColorMatrix> ds_u; // Force call will resize
  LatticeFermion X, Y;

  // Set it all up
  // Random gauge field.
  for(int mu=0; mu < Nd; mu++) { 
    gaussian(u[mu]);
    reunit(u[mu]);
  }
  
  // Fill spinors with random junk
  gaussian(X);
  gaussian(Y);


  int iter = 100;

  QDPIO::cout << "Calling Derivative : " << iter << " times " << std::endl;
  StopWatch swatch;

  swatch.reset();
  swatch.start();
  for(int i=0; i < iter; i++) { 
    deriv(u,ds_u, X, Y,0);
  }
  swatch.stop();
  QDPIO::cout << "Done in " << swatch.getTimeInSeconds() << " sec" << std::endl;

  // Possibly shutdown the machine
  QDP_finalize();

  exit(0);
}
