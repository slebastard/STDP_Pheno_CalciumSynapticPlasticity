/*
 *  stdp_structpl_names.h
 *
 *  This file is part of STDPStructplast.
 *
 *  Copyright (C) 2017, see AUTHORS.
 *
 *  STDPStructplast is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  STDPStructplast is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with STDPCalcium.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef STDPCALCIUM_NAMES_H
#define STDPCALCIUM_NAMES_H

#include "name.h"

namespace stdpcalcium
{
  /**
   * This namespace contains global Name objects. These can be used in
   * Node::get_status and Node::set_status to make data exchange
   * more efficient and consistent. Creating a Name from a std::string
   * is in O(log n), for n the number of Names already created. Using
   * predefined names should make data exchange much more efficient.
   */
  namespace names
  {
    // names specific to calcium-based plasticity model stdp_calcium_synapse
    extern const Name weight;
    extern const Name rho;
    extern const Name ca;
    extern const Name tau_ca;
    extern const Name C_pre;
    extern const Name C_post;
    extern const Name theta_p;
    extern const Name theta_d;
    extern const Name gamma_p;
    extern const Name gamma_d;
    extern const Name tau_rho;
    extern const Name S_attr;
    extern const Name sigma;
    extern const Name rho_max;
    extern const Name t_lastspike;
    extern const Name p_active; 
    extern const Name d_active; 
  }
}
#endif
