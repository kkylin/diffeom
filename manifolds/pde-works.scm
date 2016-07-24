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

;;; This stuff works!

(load "load-pde")


;;; Try something a *little* different:

(write-line '(constructing domain...))

(define make-domain
  (show-time
   (lambda ()
     (pde:make-domain-without-overlaps
      disc make-vertices make-no-extra-nodes planar-triangulate
      '(rectangular 40 15) '(spherical 20 30)))))


;;; Construct elements, as usual (cheat on the integration):

(write-line '(constructing elements...))

(show-time
 (lambda ()
   (make-domain
    imbedded-poly-laplacian
    make-triangular-imbedded-integrator
    pde:make-imbedded-poly-basis-function)))


;;; Forming matrix:

(write-line '(forming matrix...))

(define mat
  (show-time
   (lambda ()
     (combine-equations-without-overlap disc 0-function test-function))))


;;; Solve the equations:

(write-line '(relax!))

(define v
  (show-time
   (lambda ()
     (sor mat 10000 1.9))))


;;; Get a rough picture of what this looks like:

(write-line '(getting a picture of the relative error...))

(define relative-error-picture
  (show-time
   (lambda ()
     (manifold->grid 15 15 disc test-function v relative-error))))

(write-line '(getting a rough picture of the solution...))

(define solution-picture
  (show-time
   (lambda ()
     (manifold->grid 15 15 disc test-function v (lambda (val ref) val)))))


;;; We'll need to re-run these tests and save the numbers.  For now, just a
;;; brief indication (so we know what to write).

;;; What works:

;;; 1. PDE:MAKE-DOMAIN-WITH-SMALL-OVERLAPS + COMBINE-EQUATIONS-WITHOUT-OVERLAP!

;;; 2. PDE:MAKE-DOMAIN-WITHOUT-OVERLAPS, with sufficiently many nodes.

;;; And what doesn't:

;;; 1. PDE:MAKE-DOMAIN-WITH-OVERLAPS requires generating constraints, and
;;;    requires using the normal equations (which is very slow in converging).

;;; 2. PDE:MAKE-DOMAIN-WITH-LARGER-OVERLAPS generates much larger systems of
;;;    equations.  Again, using normal equations may be a bad idea.

;;; 3. PDE:MAKE-SIMPLE-DOMAIN doesn't do much better.
