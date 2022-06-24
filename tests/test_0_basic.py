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

import os
import sys
import pprint
import pytest

def test_print_cwd():
  path = os.getcwd()
  pprint.pprint(path)

def test_print_path():
  pprint.pprint(sys.path)

def test_import():
    print("import cimod")
    import cimod
    pprint.pprint(dir(cimod), compact=True)
   
def test_no_self_loop_allowed():
  import traceback
  import cimod
  import dimod
  try:
    print('BinaryQuadraticModel({}, {(2,2):5}, "SPIN")')
    cimod.BinaryQuadraticModel({}, {(2,2):5}, "SPIN")
  except RuntimeError as runtimeerror:
    pprint.pprint(runtimeerror)
    pprint.pprint(traceback.format_exc())
  else:
    print('finish (no error)')
  try:
    print('BinaryQuadraticModel({}, {(2,2):5}, cimod.SPIN)')
    cimod.BinaryQuadraticModel({}, {(2,2):5}, cimod.SPIN)
  except RuntimeError as runtimeerror:
    pprint.pprint(runtimeerror)
    pprint.pprint(traceback.format_exc())
  else:
    print('finish (no error)')
  try:
    print('BinaryQuadraticModel({}, {(2,2):5}, dimod.SPIN)')
    cimod.BinaryQuadraticModel({}, {(2,2):5}, dimod.SPIN)
  except RuntimeError as runtimeerror:
    pprint.pprint(runtimeerror)
    pprint.pprint(traceback.format_exc())
  else:
    print('finish (no error)')
  try:
    print('BinaryQuadraticModel({}, {(2,2):5}, cimod.cxxcimod.Vartype.SPIN)')
    cimod.BinaryQuadraticModel({}, {(2,2):5}, cimod.cxxcimod.Vartype.SPIN)
  except RuntimeError as runtimeerror:
    pprint.pprint(runtimeerror)
    pprint.pprint(traceback.format_exc())
  else:
    print('finish (no error)')
