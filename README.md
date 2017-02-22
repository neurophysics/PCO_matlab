# Phase Coupling Optimization

This implements Phase Coupling Optimization in Matlab

A citable reference is not ready yet, but on the way.


***************************************************************************
Disclosure:                                                       
-----------                                                       
This software comes as it is - there might be errors at runtime and results
might be wrong although the code was tested and did work as expected. Since
resuluts might be wrong you must absolutely not use this software for a
medical purpuse - decisions concerning diagnosis, treatment or prophylaxis
***************************************************************************

Usage:
------
[vlen, wy] = PCOa(a, y)
[vlen, wy] = PCOa(a, y, nu)
[vlen, wy] = PCOa(a, y, nu, bestof)

Standard Parameters:
------------------
nu = 1
bestof = 15

Inputs:
-------
a - (column vector) amplitudes
y - (2d array, complex) analytic representation of signal,
    channels x datapoints
num - (int > 0) - determines the number of filters that will be
                  derived. This depends also on the rank of Y,
                  the final number of filters will be min
                  ([num, rank(y)]), defaults to 1
bestof (int > 0) - the number of restarts for the optimization of the
                   individual filters. The best filter over all
                   these restarts with random initializations is
                   chosen, defaults to 15.

Outputs:
-------
vlen - row vector - the length of the mean vector for each filter
wy - 2d array - the filters for Y, each filter is in a column of Wy
***************************************************************************

License:
--------
Copyright (c) 2017 Gunnar Waterstraat & Iara de Almeida Ivo

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Author & Contact
----------------
Written by Gunnar Waterstraat & Iara de Almeida Ivo
email: gunnar[dot]waterstraat[at]charite.de

