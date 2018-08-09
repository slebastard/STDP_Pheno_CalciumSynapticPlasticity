#ifndef STDP_CONNECTION_CALCIUM
#define STDP_CONNECTION_CALCIUM

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

namespace nest
{
  class STDPCalcium : public StaticConnection
  {
    public:
      STDPCalcium () {}
      STDPCalcium (const STDPCalcium &) {}
      ~STDPCalcium () {}

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
      void send (Event & e, double_t t_lastspike, const CommonSynapseProperties & cp);


      void set_weight( double w )
      {
        weight_ = w;
        // Use inverse transfer to set conresponding rho
      }

      void set_calcium( double c )
      {
        ca_ = c;
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

      double transfer_( double rho )
      {
        return erf((S_attr_ - rho)/sigma_);
      }

      double invtransfer_( double w )
      {
        return std::min(std::max(0, S_attr_ - 2*std::pow(sigma_,2.0)*erfc_inv(w)), rho_max_);
      }

      double calcium_update_prepost_( double dt )
      {
        return ca_*std::exp(-dt/tau_ca_) + C_post_;
      }

      double calcium_update_postpre_( double dt )
      {
        return ca_*std::exp(-dt/tau_ca_) + C_pre_;
      }

      double update_rho_( double time_pot, double time_dep )
      {
        return ;
      }
  };

  inline void STDPCalcium::send (Event & e, double_t t_lastspike, const CommonSynapseProperties &)
  {
    update_dynamics();

    e.set_receiver(*target_);
    e.set_weight(weight_);
    e.set_delay(delay_);
    e.set_rport(rport_);
    e();
  }

template < typename targetidentifierT >
inline void
STDPConnection< targetidentifierT >::send( Event& e,
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
    ca_ = calcium_update_prepost_( ca_delay );

    t_max_p_duration = p_active_ * tau_ca_*std::log(ca_/theta_p_);
    t_close_p = post_previousspike + std::min(p_active_ * (post_spike - post_previousspike), t_max_p_duration);

    t_max_d_duration = d_active_ * tau_ca_*std::log(ca_/theta_d_);
    t_close_d = post_previousspike + std::min(d_active_ * (post_spike - post_previousspike), t_max_d_duration);

    p_active_ = (ca_ >= theta_p_);
    d_active_ = (ca_ >= theta_d_);

    rho_ = update_rho_( t_close_p, t_close_d );
    weight_ = transfer_( rho_ );
  }


  target->get_history( t_lastspike_ - dendritic_delay,
    t_spike - dendritic_delay,
    &start,
    &finish );

  ca_ = calcium_update_postpre_( t_spike - post_lastspike );
  p_active_ = (ca_ >= theta_p_);
  d_active_ = (ca_ >= theta_d_);
  while ( finish != start )
  {
    post_spike = finish->t_;
    dt_postpre = post_spike + dendritic_delay - t_spike;
    --finish;
    // Postpre sweep
    assert(dt_postpre<=0);

    rho_ = update_rho_( dt_postpre );
    weight_ = transfer_( rho_ );
  }

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
STDPCalcium< targetidentifierT >::STDPConnection()
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
STDPCalcium< targetidentifierT >::STDPConnection(
  const STDPConnection< targetidentifierT >& rhs )
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
STDPCalcium< targetidentifierT >::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< double >( d, names::rho_, rho_ );
  def< double >( d, names::ca_, ca_ );
  def< double >( d, names::tau_ca_, tau_ca_ );
  def< double >( d, names::C_pre_, C_pre_ );
  def< double >( d, names::C_post_, C_post_ );
  def< double >( d, names::theta_p_, theta_p_ );
  def< double >( d, names::theta_d_, theta_d_ );
  def< double >( d, names::gamma_p_, gamma_p_ );
  def< double >( d, names::gamma_d_, gamma_d_ );
  def< double >( d, names::tau_rho_, tau_rho_ );
  def< double >( d, names::S_attr_, S_attr_ );
  def< double >( d, names::sigma_, sigma_ );
  def< double >( d, names::rho_max_, rho_max_ );
  def< double >( d, names::t_lastspike_, t_lastspike_ );
  def< bool >( d, names::p_active_, p_active_ ); 
  def< bool >( d, names::d_active_, d_active_ ); 
  def< long >( d, names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
STDPCalcium< targetidentifierT >::set_status( const DictionaryDatum& d,
  ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
  updateValue< double >( d, names::rho_, rho_ );
  updateValue< double >( d, names::ca_, ca_ );
  updateValue< double >( d, names::tau_ca_, tau_ca_ );
  updateValue< double >( d, names::C_pre_, C_pre_ );
  updateValue< double >( d, names::C_post_, C_post_ );
  updateValue< double >( d, names::theta_p_, theta_p_ );
  updateValue< double >( d, names::theta_d_, theta_d_ );
  updateValue< double >( d, names::gamma_p_, gamma_p_ );
  updateValue< double >( d, names::gamma_d_, gamma_d_ );
  updateValue< double >( d, names::tau_rho_, tau_rho_ );
  updateValue< double >( d, names::S_attr_, S_attr_ );
  updateValue< double >( d, names::sigma_, sigma_ );
  updateValue< double >( d, names::rho_max_, rho_max_ );
  updateValue< double >( d, names::t_lastspike_, t_lastspike_ );
  updateValue< bool >( d, names::p_active_, p_active_ ); 
  updateValue< bool >( d, names::d_active_, d_active_ ); 

  // check if weight_ and Wmax_ has the same sign
  if ( not( ( ( rho_ >= 0 ) - ( rho_ < 0 ) )
         == ( ( rho_max_ >= 0 ) - ( rho_max_ < 0 ) ) ) )
  {
    throw BadProperty( "rho and rho_max must have same sign." );
  }
}


} // namespace nest

#endif /* #ifndef STDP_CONNECTION_CALCIUM */