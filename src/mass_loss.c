/**
 * @file    mass_loss.c
 * @brief   
 * @author  
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Dan Tamayo, Hanno Rein
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Mass Modifications$     // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script). 
 * 
 * ======================= ===============================================
 * Authors                 
 * Implementation Paper    *In progress*
 * Based on                None
 * C Example               :ref:`c_example_modify_mass`
 * Python Example          `ModifyMass.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/ModifyMass.ipynb>`_.
 * ======================= ===============================================
 * 
 * This adds exponential mass growth/loss to individual particles every timestep.
 * Set particles' ``tau_mass`` parameter to a negative value for mass loss, positive for mass growth.
 * 
 * **Effect Parameters**
 * 
 * *None*
 * 
 * **Particle Parameters**
 * 
 * Only particles with their ``tau_mass`` parameter set will have their masses affected.
 * 
 * ============================ =========== =======================================================
 * Name (C type)                Required    Description
 * ============================ =========== =======================================================
 * tau_mass (double)            Yes         e-folding mass loss (<0) or growth (>0) timescale    
 * ============================ =========== =======================================================
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void rebx_mass_loss(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){