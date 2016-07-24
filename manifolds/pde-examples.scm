#|

This file is part of DIFFEOM ("DIFFerential Equations
On Manifolds"), a system for solving differential
equations on manifolds.

Copyright (C) 2016 by Kevin K Lin
<klin@math.arizona.edu>

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

;;; We won't get any farther on this endeavor without defining some good
;;; examples.  A good example to use first is Laplace's equation, since the
;;; answers are easy to check.  Another option is to define something on the
;;; torus or the sphere -- How does one define the Laplacian in that case?
;;; Laplacian = (d + d*)^2 = d*d + dd*.  On 0-forms (smooth functions) d* = 0,
;;; so it's just d*d.  How does one compute the adjoint locally?


;;; Laplace's equation in a compact region of the plane (with boundary),
;;; covered by more than one coordinate system.

;;; First, let the domain be the unit closed disc, with spherical coordinates
;;; on the boundary:

(define disc
  (make-ball 2 make-spherical-sphere))


;;; Define some Laplacian operators (on different basis functions):

(define poly-disc-laplacian
  (make-operator
   disc
   (operator:pull-back-poly-op
    poly-gradient
    poly-gradient
    (lambda (v w) (basis:scalar* -1 (basis:dot v w))))))

(define imbedded-poly-laplacian
  (make-operator
   disc
   (operator:imbedded-poly-op
    poly-gradient
    poly-gradient
    (lambda (v w) (basis:scalar* -1 (basis:dot v w))))))

(define real-disc-laplacian
  (make-operator
   disc
   (operator:pull-back-real-op
    real-gradient
    real-gradient
    (lambda (v w) (basis:scalar* -1 (basis:dot v w))))))


;;; Some other things that are useful as test solutions:

(define x-coord-map
  (compose vector-first node:get-point))

(define y-coord-map
  (compose vector-second node:get-point))

(define (test-function node)
  (let ((x (x-coord-map node))
	(y (y-coord-map node)))
    (- (square x) (square y))))


;;; Useful definitions to have (for debugging purposes):

(define atlas (manifold:get-finite-atlas disc))
(define c1 (car atlas))
(define c2 (cadr atlas))
(define c3 (caddr atlas))


;;; The null equation corresponds to the node in row 228 of MAT.  It is the
;;; node at the center of the disc.  Why are its coefficients 0?  We should
;;; never get null equations because the diagonal terms should at least be
;;; non-zero.

;;; It may appear that the problem lies with the fact that differential
;;; operator operates on (inexact, interpolated) pull-backs of basis functions
;;; back onto the imbedded-manifold.

;;; But why should this cause problems on C3, which is really a subset of R^2?
;;; What really must be happening is that Galerkin's method doesn't work unless
;;; one performs the integration by parts, which is not exactly kosher because
;;; the basis functions are only piecewise-differentiable.  So a variational
;;; principle must be really what's at work in finite elements...

;;; Using first-order operators and inner products seems to get rid of the
;;; problem.  This is because while the basis functions are C^2 over elements,
;;; they are only C^0 across edges.  Thus, differential operator we can safely
;;; apply to basis functions have at most order 1, and it is necessary to split
;;; the operator using integration by parts.

;;; Most of the error is probably coming from truncation errors in computing
;;; the differential operator, since we cut off higher-order terms when
;;; applying coordinate transformations.  The other possible source of error is
;;; the constraint.  First thing to try, then, is to implement a nice
;;; multidimensional numerical integrator.
