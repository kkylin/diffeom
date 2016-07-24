#|

This file is part of DIFFEOM, a system for solving
differential equations on manifolds.

Copyright (C) 2016 by Kevin K Lin
<kkylin@alum.mit.edu>

This program is free software; you can redistribute
it and/or modify it under the terms of the GNU
General Public License as published by the Free
Software Foundation; either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General
Public License along with this program; if not, write
to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301 USA.

|#

(load "util")            ;; Useful helper routines.
(load "util-too")        ;; Useful routines compatible with ScmUtils.
(load "matlib")          ;; Everyone needs the matrix library.
(load "sparse")          ;; Sparse matrices.
(load "relax")           ;; Relaxation.
(load "dyntable")        ;; Dynamic tables.
(load "debug")           ;; Graphics & stuff.
(load "delaunay")        ;; Delaunay triangulation.
(load "delaux")          ;; Auxiliary routines for Delaunay.
(load "edge")            ;; Edge algebra junk for Delaunay.
(load "fem")             ;; The main finite-element code.
(load "nodes")           ;; Definition of nodes.
(load "basis")           ;; Basis functions.
(load "2d-poly-basis")   ;; Polynomial basis functions in two variables.
(load "2d-real-diff")    ;; Real differential operators on real functions.
(load "2d-real-basis")   ;; Polynomial basis functions in two variables.
(load "2d-trapezoid")    ;; Numerical integration using trapezoidal rule.
(load "operators")       ;; Tools for differential operators.
(load "2d-operators")    ;; Examples.
(load "2d-domains")      ;; Making domains and boundaries, etc.
(load "2d-examples")     ;; Some examples...
