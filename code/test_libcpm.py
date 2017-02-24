# -*-coding:Utf-8 -*
# ====================================================================
# Copyright 2017 Clement Ranc
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied. See the License for the specific language governing
# permissions and limitations under the License.
# ====================================================================
# Standard packages
# ====================================================================
import sys
import numpy as np
# ====================================================================
# Non-standard packages
# ====================================================================
import libcpm
# ====================================================================
#   Functions
# ====================================================================
def test_magnifcalc(**kwargs):

    s = np.linspace(0,10, 100)
    x = np.linspace(0, 10, 100)
    y = np.linspace(0, 10, 100)
    list = [[0,0,0,0,0,0,0], s.tolist(), x.tolist(), y.tolist()]
    result = libcpm.linear_least_squares_wrap(list)

    return result


# ====================================================================
# Main function
# ====================================================================
if __name__=="__main__":
    test_magnifcalc()










