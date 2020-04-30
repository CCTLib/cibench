/*! \file
 *  \brief Plaquette gauge action as sum of characters
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/character_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"
#include "io/aniso_io.h"

namespace Chroma
{
 
  namespace CharacterGaugeActEnv 
  { 
    namespace
    {
      GaugeAction< multi1d<LatticeColorMatrix>, 
		   multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
								 const std::string& path) 
      {
	return new GaugeAct(CreateGaugeStateEnv::reader(xml, path), 
			    Params(xml, path));
      }
      
      const std::string name = "CHARACTER_GAUGEACT";

      //! Local registration flag
      static bool registered = false;
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheGaugeActFactory::Instance().registerObject(name, createGaugeAct);
	registered = true;
      }
      return success;
    }


    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      XMLReader paramtop(xml_in, path);

      try 
      {	
	read(paramtop, "beta_F", beta_F);
	read(paramtop, "beta_A", beta_A);
	read(paramtop, "beta_S", beta_S);
      }
      catch( const std::string& e ) { 
	QDPIO::cerr << "Error reading XML: " <<  e << std::endl;
	QDP_abort(1);
      }
    }


    //! Compute the action
    Double GaugeAct::S(const Handle< GaugeState<P,Q> >& state) const
    {
      // Action at the site level
      multi2d<LatticeColorMatrix> plq;
      this->siteAction(plq, state);

      // Total action
      // Fundamental part
      Double act_F  = zero;
      Double act_A  = zero;
      Double act_S  = zero;
      Double three = Nc;
      Double nine = Nc*Nc ;
      Double sext = Nc*Nc+Nc ;

      for(int mu=1; mu < Nd; ++mu)
      {
	for(int nu=0; nu < mu; ++nu)
	{
	  LatticeComplex plaq=trace(plq[mu][nu]);
	  // Sum over plaquettes
	  act_F += sum(real(plaq)-three);
	  act_A += sum(localNorm2(plaq) - nine);
	  act_S += sum(real(trace(plq[mu][nu]*plq[mu][nu])+plaq*plaq)-sext);
	}
      }

      // Normalize
      Real act = -(param.beta_F / Real(Nc)) * act_F - (param.beta_A/Real(Nc*Nc))  * act_A - 0.5*param.beta_S*(act_S/sext);

      return act;
    }
 

    //! Compute the plaquette
    void GaugeAct::siteAction(multi2d<LatticeColorMatrix>& plq, const Handle< GaugeState<P,Q> >& state) const
    {
      START_CODE();

      // Initialize
      plq.resize(Nd,Nd);
      plq = zero;

      // Handle< const GaugeState<P,Q> > u_bc(createState(u));
      // Apply boundaries
      const multi1d<LatticeColorMatrix>& u = state->getLinks();

      // Compute the average plaquettes
      for(int mu=1; mu < Nd; ++mu)
      {
	for(int nu=0; nu < mu; ++nu)
	{
	  /* tmp_0 = u(x+mu,nu)*u_dag(x+nu,mu) */
	  /* tmp_1 = tmp_0*u_dag(x,nu)=u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu) */
	  /* wplaq_tmp = tr(u(x,mu)*tmp_1=u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu)) */
	  plq[mu][nu] += (u[mu]*shift(u[nu],FORWARD,mu)*adj(shift(u[mu],FORWARD,nu))*adj(u[nu])) ; 
	  
	  // Keep a copy
	  plq[nu][mu] = plq[mu][nu];
	}
      }


      END_CODE();
    }
 

    //! Compute staple
    /*! Default version. Derived class should override this if needed. */
    void GaugeAct::staple(LatticeColorMatrix& result,
			  const Handle< GaugeState<P,Q> >& state,
			  int mu, int cb) const
    {
      QDPIO::cerr << __func__ << ": staple not possible\n";
      QDP_abort(1);
    }


    //! Compute dS/dU
    void GaugeAct::deriv(multi1d<LatticeColorMatrix>& ds_u,
			 const Handle< GaugeState<P,Q> >& state) const
    {
      START_CODE();

      LatticeColorMatrix tmp_1;
      LatticeColorMatrix tmp_2;

      const multi1d<LatticeColorMatrix>& u = state->getLinks();
      multi1d<LatticeColorMatrix> deriv_fun(Nd); 
      multi1d<LatticeColorMatrix> deriv_adj(Nd);
      multi1d<LatticeColorMatrix> deriv_sex(Nd);
      ds_u.resize(Nd);
      for(int mu=0; mu < Nd; mu++) 
      {
	deriv_fun[mu] = zero ;
	deriv_adj[mu] = zero ;
	deriv_sex[mu] = zero ;
	for(int nu = 0; nu < Nd; nu++) 
	{ 
	  if (mu == nu) continue;

	  LatticeColorMatrix tmp_1 = shift(u[nu], FORWARD, mu);
	  LatticeColorMatrix tmp_2 = shift(u[mu], FORWARD, nu);

	  LatticeColorMatrix up_plq   = u[mu]*tmp_1*adj(tmp_2)*adj(u[nu]);
	  LatticeColorMatrix down_plq = u[mu]*shift(adj(tmp_1)*adj(u[mu])*u[nu],BACKWARD,nu);

	  deriv_fun[mu] += up_plq + down_plq ;
	  
	  deriv_adj[mu] += up_plq*conj(trace(up_plq)) + 
	                   down_plq*conj(trace(down_plq)) ;

	  deriv_sex[mu] += up_plq*up_plq + down_plq*down_plq +
	                   up_plq*trace(up_plq) + 
	                   down_plq*trace(down_plq) ;
	
	}// nu
	// Fold in the normalization from the action
	ds_u[mu]  = (-param.beta_F/Real(2*Nc)        ) * deriv_fun[mu];
	ds_u[mu] += (-param.beta_A/Real(Nc*Nc)       ) * deriv_adj[mu];
	ds_u[mu] += (-param.beta_S/Real(2*(Nc*Nc+Nc))) * deriv_sex[mu];
      }// mu

      // Zero the force on any fixed boundaries
      getGaugeBC().zero(ds_u);

      END_CODE();
    }


  }//CharacterGaugeActEnv

} // Chroma
