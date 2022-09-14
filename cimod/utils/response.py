# Copyright 2022 Jij Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import annotations

import dimod

def get_state_and_energy(
    model, result_state, offset=0, model_variables=[]
) -> tuple[dict, float]:
    """get converted state and energy.
    This function receives raw array of spins or binaries.
    If vartype of model and the vartype of the raw array are different, the raw array is automatically converted to the vartype of model with any offset shift.

    Args:
        model: cimod model (BinaryQuadraticModel or BinaryPolynomialModel)
        result_state (list): states of spins or binaries
        offset (float): offset added to returned energy
        model_variables (Optional[list]): list of variables

    Returns:
        tuple[dict, float]: labeled states and corresponding energy

    Examples:


    """
    temp_state = {}

    if len(model_variables) != 0:
        variables = model_variables
    else:
        variables = model.variables

    def _convert_data(d):
        if d == 0 and model.vartype == dimod.SPIN:
            return -1
        elif d == -1 and model.vartype == dimod.BINARY:
            return 0
        else:
            return d

    for num in range(len(model.variables)):
        temp_state[variables[num]] = _convert_data(result_state[num])

    return temp_state, model.energy(temp_state) + offset
