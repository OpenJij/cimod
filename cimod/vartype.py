# Copyright 2020 Jij Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import cxxcimod
import dimod

SPIN = dimod.SPIN
BINARY = dimod.BINARY
Vartype = dimod.Vartype

def to_cxxcimod(var_type):
    # convert to cxxcimod type
    if isinstance(var_type, cxxcimod.Vartype):
        return var_type
    
    var_type = dimod.as_vartype(var_type)
    if var_type == dimod.SPIN:
        return cxxcimod.Vartype.SPIN
    if var_type == dimod.BINARY:
        return cxxcimod.Vartype.BINARY
