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

def disabled(func):
    def wrapper(*args, **kwargs):
        raise NotImplementedError("The function {} is disabled.".format(func.__name__))

    return wrapper


def recalc(func):
    def wrapper(self, *args, **kwargs):
        self._re_calculate = True
        self._re_calculate_indices = True
        return func(self, *args, **kwargs)

    return wrapper
