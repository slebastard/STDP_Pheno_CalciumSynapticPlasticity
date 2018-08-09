/*
 *  stdp_connection_calcium.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef STDP_CONNECTION_CALCIUM
#define STDP_CONNECTION_CALCIUM

/* BeginDocumentation
  Name: stdp_connection_calcium - Synapse type for spike-timing dependent
   plasticity, based on calcium-history, with stability transfer function
   (individual synapse is assumed to be bistable)

  Description:
   stdp_connection_calcium is a connector to create synapses with a spike time
   dependent plasticity based on a calcium history.
   The model is a modified version of [1].
   It integrates a phenomenological dynamics for a phosphorylation variable, on which the synaptic weight directly depends
   (see Graupner & Brunel 2007, http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.0030221)
   
   The model relies on the following equations:
   tau_ca * d[ca]/dt = -ca + sum[C_pre*dirac[t_pre] + c_post*dirac[t_post]]
   tau_rho * d[rho]/dt = gamma_p*(rho_max - rho)*(ca > theta_p) - gamma_d*rho*(ca > theta_d)
   tau_w * d[weight]/dt = transfer(rho) - weight

  Variables:
   weight    double    The synaptic weight, between 0 and 1
   rho       double    A phosphorylation proxy variable that directly determines the weight
   ca        double    The current calcium level at the synapse

  Parameters:
   tau_ca    double    The time constant of calcium decay at the synapse
   C_pre     double    The amount of calcium added by a presynaptic spike
   C_post    double    The amount of calcium added by a postsynaptic spike
   theta_p   double    The threshold that the calcium level must exceed for potentiation to take place
   theta_d   double    The threshold that the calcium level must exceed for depression to take place
   gamma_p   double    Controls the speed of potentiation of the phosphorylation variable
   gamma_d   double    Controls the speed of depression of the phosphorylation variable
   tau_rho   double    Time constant for phosphorylation
   S_attr    double    Phosphorylation level for the saddle node bifurcation that determines bistability of the phosphorylation var at basal calcium level
   sigma     double    Noise level
   rho_max   double    Maximum phosphorylation level
   p_active_ bool      Active potentiation flag (ca > theta_p)?
   d_active_ bool      Active depression flag (ca > theta_d)?

  Transmits: SpikeEvent

  References:
   [1]  Graupner M, Brunel N (2012)
        Calcium-based plasticity model explains sensitivity of synaptic changes to spike pattern, rate, and dendritic location.
        Proc Natl Acad Sci U S A. 2012 Mar 6; 109(10): 3991–3996.
        https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3309784/

   [2]  Graupner M, Brunel N (2007)
        STDP in a Bistable Synapse Model Based on CaMKII and Associated Signaling Pathways.
        PLoS Comput Biol 3(11): e221.
        https://doi.org/10.1371/journal.pcbi.0030221

  FirstVersion: August 2018
  Author: Simon Lebastard
*/

// C++ includes:
#include <cmath>

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connection.h"
#include "connector_model.h"
#include "event.h"

// Includes from sli:
#include "dictdatum.h"
#include "dictutils.h"

// Including boost for inverse erfc function
#include <boost/math/special_functions/erf.hpp>

namespace stdpcalcium
{
  class STDPCalciumConnector : public StaticConnection
  {
    public:
      STDPCalciumConnector () {}
      STDPCalciumConnector (const STDPCalciumConnector &) {}
      ~STDPCalciumConnector () {}

      /**
       * Get all properties of this connection and put them into a dictionary.
       */
      void get_status( DictionaryDatum& d ) const;

      /**
       * Set properties of this connection from the values given in dictionary.
       */
      void set_status( const DictionaryDatum& d, ConnectorModel& cm );

      /**
       * Send an event to the receiver of this connection.
       * \param e The event to send
       * \param cp common properties of all synapses (empty).
       */
      void send (Event& e, thread t, const CommonSynapseProperties& cp);

      void
      set_weight( double w )
      {
        weight_ = w;
      }

    private:
      double weight_;
      double rho_;
      double ca_;

      const double tau_ca_;
      const double C_pre_;
      const double C_post_;
      const double theta_p_;
      const double theta_d_;
      const double gamma_p_;
      const double gamma_d_;
      const double tau_rho_;
      const double S_attr_;
      const double sigma_;
      const double rho_max_;

      double t_lastspike_;
      bool p_active_;
      bool d_active_;

      double transfer_( double rho );

      double invtransfer_( double w );

      void calcium_update_prepost_( double dt );

      void calcium_update_postpre_( double dt );

      void update_rho_( double time_pot, double time_dep );
  };

  inline double transfer_( double rho )
  {
    return std::erfc((S_attr_ - rho)/sigma_);
  }

  inline double invtransfer_( double w )
  {
    return std::min(std::max(0, S_attr_ - 2*std::pow(sigma_,2.0)*erfc_inv(w)), rho_max_);
  }
  
  inline void calcium_update_prepost_( double dt )
  {
    ca_ = ca_*std::exp(-dt/tau_ca_) + C_post_;
  }
  
  inline void calcium_update_postpre_( double dt )
  {
    ca_ = ca_*std::exp(-dt/tau_ca_) + C_pre_;
  }
  
  inline void update_rho_( double time_pot, double time_dep )
  {
    double t_pure_d = t_close_d - t_close_p;
    assert( t_pure_d >= 0 );
  
    rho_ = ((rho_ - rho_max_*gamma_p_/(gamma_p_+gamma_d_))*std::exp(-(t_close_p)*(gamma_p_+gamma_d_)/tau_rho_) + rho_max_*gamma_p_/(gamma_p_+gamma_d_));
    rho_ = rho_.* std::exp(-gamma_d_*(t_pure_d)/tau_rho_);
  
    weight_ = transfer_( rho_ );
  }

  template < typename targetidentifierT >
  inline void
  STDPCalciumConnector< targetidentifierT >::send( Event& e,
    thread t,
    const CommonSynapseProperties& )
  {
    // synapse STDP depressing/facilitation dynamics
    double t_spike = e.get_stamp().get_ms();
  
    // use accessor functions (inherited from Connection< >) to obtain delay and
    // target
    Node* target = get_target( t );
    double dendritic_delay = get_delay();
  
    // get spike history in relevant range (t1, t2] from post-synaptic neuron
    std::deque< histentry >::iterator start;
    std::deque< histentry >::iterator finish;
  
    // For a new synapse, t_lastspike_ contains the point in time of the last
    // spike. So we initially read the
    // history(t_last_spike - dendritic_delay, ..., T_spike-dendritic_delay]
    // which increases the access counter for these entries.
    // At registration, all entries' access counters of
    // history[0, ..., t_last_spike - dendritic_delay] have been
    // incremented by Archiving_Node::register_stdp_connection(). See bug #218 for
    // details.
    target->get_history( t_lastspike_ - dendritic_delay,
      t_spike - dendritic_delay,
      &start,
      &finish );
    // facilitation due to post-synaptic spikes since last pre-synaptic spike
    double post_spike;
    double post_previousspike;
    double post_lastspike;
    double dt_prepost;
    double dt_postpre;
    double ca_delay;

    double t_max_p_duration;
    double t_close_p;
    double t_max_d_duration;
    double t_close_d;

    post_lastspike = finish->t_;
    post_previousspike = t_lastspike_;
  
    while ( start != finish )
    {
      post_spike = start->t_;
      dt_prepost = post_spike + dendritic_delay - t_lastspike_;
      ++start;
      ca_delay = post_spike - post_previousspike;
      post_previousspike = post_spike;
  
      // Prepost sweep
      assert(dt_prepost>=0);
  
      t_max_p_duration = p_active_ * tau_ca_*std::log(ca_/theta_p_);
      t_close_p = std::min(p_active_ * (post_spike - post_previousspike), t_max_p_duration);
  
      t_max_d_duration = d_active_ * tau_ca_*std::log(ca_/theta_d_);
      t_close_d = std::min(d_active_ * (post_spike - post_previousspike), t_max_d_duration);
  
      calcium_update_prepost_( ca_delay );
      p_active_ = (ca_ >= theta_p_);
      d_active_ = (ca_ >= theta_d_);
  
      update_rho_( t_close_p, t_close_d );
    }
  
    t_max_p_duration = p_active_ * tau_ca_*std::log(ca_/theta_p_);
    t_close_p = std::min(p_active_ * (t_spike - post_previousspike), t_max_p_duration);
  
    t_max_d_duration = d_active_ * tau_ca_*std::log(ca_/theta_d_);
    t_close_d = std::min(d_active_ * (t_spike - post_previousspike), t_max_d_duration);
  
    ca_ = calcium_update_postpre_( t_spike - post_lastspike );
    p_active_ = (ca_ >= theta_p_);
    d_active_ = (ca_ >= theta_d_);
    rho_ = update_rho_( t_close_p, t_close_d );
  
  
    e.set_receiver( *target );
    e.set_weight( weight_ );
    // use accessor functions (inherited from Connection< >) to obtain delay in
    // steps and rport
    e.set_delay( get_delay_steps() );
    e.set_rport( get_rport() );
    e();
  
    t_lastspike_ = t_spike;
  }
  
  
  
  template < typename targetidentifierT >
  STDPCalciumConnector< targetidentifierT >::STDPCalciumConnector()
    : ConnectionBase()
    , weight_(0.0)
    , rho_(0.0)
    , ca_(0.0)
    , tau_ca_(80.0)
    , C_pre_(0.40)
    , C_post_(0.84)
    , theta_p_(1.08)
    , theta_d_(1.00)
    , gamma_p_(120.0)
    , gamma_d_(200.0)
    , tau_rho_(100000.0)
    , S_attr_(40.0)
    , sigma_(25.0)
    , rho_max_(200.0)
    , t_lastspike_( 0.0 )
    , p_active_( false )
    , d_active_( false )
  {
  }
  
  template < typename targetidentifierT >
  STDPCalciumConnector< targetidentifierT >::STDPCalciumConnector(
    const STDPCalciumConnector< targetidentifierT >& rhs )
    : ConnectionBase( rhs )
    , weight_( rhs.weight_ )
    , rho_( rhs.rho_ )
    , ca_( rhs.ca_ )
    , tau_ca_( rhs.tau_ca_ )
    , C_pre_( rhs.C_pre_ )
    , C_post_( rhs.C_post_ )
    , theta_p_( rhs.theta_p_ )
    , theta_d_( rhs.theta_d_ )
    , gamma_p_( rhs.gamma_p_ )
    , gamma_d_( rhs.gamma_d_ )
    , tau_rho_( rhs.tau_rho_ )
    , S_attr_( rhs.S_attr_ )
    , sigma_( rhs.sigma_ )
    , rho_max_( rhs.rho_max_ )
    , t_lastspike_( rhs.t_lastspike_ )
    , p_active_( rhs.p_active_ )
    , d_active_( rhs.d_active_ )
  {
  }
  
  template < typename targetidentifierT >
  void
  STDPCalciumConnector< targetidentifierT >::get_status( DictionaryDatum& d ) const
  {
    ConnectionBase::get_status( d );
    def< double >( d, names::weight, weight_ );
    def< double >( d, names::rho, rho_ );
    def< double >( d, names::ca, ca_ );
    def< double >( d, names::tau_ca, tau_ca_ );
    def< double >( d, names::C_pre, C_pre_ );
    def< double >( d, names::C_post, C_post_ );
    def< double >( d, names::theta_p, theta_p_ );
    def< double >( d, names::theta_d, theta_d_ );
    def< double >( d, names::gamma_p, gamma_p_ );
    def< double >( d, names::gamma_d, gamma_d_ );
    def< double >( d, names::tau_rho, tau_rho_ );
    def< double >( d, names::S_attr, S_attr_ );
    def< double >( d, names::sigma, sigma_ );
    def< double >( d, names::rho_max, rho_max_ );
    def< double >( d, names::t_lastspike, t_lastspike_ );
    def< bool >( d, names::p_active, p_active_ ); 
    def< bool >( d, names::d_active, d_active_ ); 
    def< long >( d, names::size_of, sizeof( *this ) );
  }
  
  template < typename targetidentifierT >
  void
  STDPCalciumConnector< targetidentifierT >::set_status( const DictionaryDatum& d,
    ConnectorModel& cm )
  {
    ConnectionBase::set_status( d, cm );
    updateValue< double >( d, names::weight, weight_ );
    updateValue< double >( d, names::rho, rho_ );
    updateValue< double >( d, names::ca, ca_ );
    updateValue< double >( d, names::tau_ca, tau_ca_ );
    updateValue< double >( d, names::C_pre, C_pre_ );
    updateValue< double >( d, names::C_post, C_post_ );
    updateValue< double >( d, names::theta_p, theta_p_ );
    updateValue< double >( d, names::theta_d, theta_d_ );
    updateValue< double >( d, names::gamma_p, gamma_p_ );
    updateValue< double >( d, names::gamma_d, gamma_d_ );
    updateValue< double >( d, names::tau_rho, tau_rho_ );
    updateValue< double >( d, names::S_attr, S_attr_ );
    updateValue< double >( d, names::sigma, sigma_ );
    updateValue< double >( d, names::rho_max, rho_max_ );
    updateValue< double >( d, names::t_lastspike, t_lastspike_ );
    updateValue< bool >( d, names::p_active, p_active_ ); 
    updateValue< bool >( d, names::d_active, d_active_ ); 
  
    // check if weight_ and Wmax_ has the same sign
    if ( not( ( ( rho_ >= 0 ) - ( rho_ < 0 ) )
           == ( ( rho_max_ >= 0 ) - ( rho_max_ < 0 ) ) ) )
    {
      throw BadProperty( "rho and rho_max must have same sign." );
    }
  }

} // namespace nest

#endif /* #ifndef STDP_CONNECTION_CALCIUM */