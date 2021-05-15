from typing import Tuple
import cimod

def get_state_and_energy(model, result_state, offset=0, model_variables=[]) -> Tuple[dict, float]:
    """get converted state and energy.
    This function receives raw array of spins or binaries.
    If vartype of model and the vartype of the raw array are different, the raw array is automatically converted to the vartype of model with any offset shift.

    Args:
        model: cimod model (BinaryQuadraticModel or BinaryPolynomialModel)
        result_state (list): states of spins or binaries
        offset (float): offset added to returned energy
        model_variables (Optional[list]): list of variables

    Returns:
        Tuple[dict, float]: labeled states and corresponding energy

    Examples:
        
        
    """
    temp_state = {}

    if len(model_variables) != 0:
        variables = model_variables
    else:
        variables = model.variables

    def _convert_data(d):
        if d == 0 and model.vartype == cimod.SPIN:
            return -1
        elif d == -1 and model.vartype == cimod.BINARY:
            return 0
        else:
            return d

    for num in range(len(model.variables)):
        temp_state[variables[num]] = _convert_data(result_state[num])

    return temp_state, model.energy(temp_state)+offset

