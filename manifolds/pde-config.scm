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

;;; This file defines procedures that need to be easily modifiable:


;;; Methods for generating equations from triangulated domains:

(define combine-equations-without-overlap
  (pde:equation-maker merge-equations))

(define combine-equations-with-overlap
  (pde:equation-maker
   (append-constraint-equations make-ordered-boundary-constraints)))

(define combine-equations-with-overlap1
  (pde:equation-maker
   (append-constraint-equations make-all-constraints)))

(define combine-equations-with-overlap2
  (pde:equation-maker
   (append-constraint-equations make-all-ordered-constraints)))

(define combine-equations-using-CMPGRD
  (pde:equation-maker cmpgrd:combine-equations))


;;; Methods for triangulating domains:

(define pde:make-domain-without-overlaps
  (pde:domain-maker generate-node-lists exact-overlap))

(define pde:make-domain-with-small-overlaps
  (pde:domain-maker copy-between-node-lists exact-overlap))

(define pde:make-domain-with-overlaps
  (pde:domain-maker make-nodes-for-each-chart reduce-overlap))

(define pde:make-domain-with-larger-overlaps
  (pde:domain-maker make-nodes-for-each-chart extended-overlap))

(define pde:make-simple-domain
  (pde:domain-maker make-nodes-for-each-chart do-nothing-to-complex))
