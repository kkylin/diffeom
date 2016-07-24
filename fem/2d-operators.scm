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

;;; Some examples:

(declare (usual-integrations))


;;; The two-dimensional Laplcian:

(define (laplacian nodes)
  (make-operator
   poly-gradient
   poly-gradient
   (lambda (v w) (basis:scalar* -1 (basis:dot v w)))))


;;; The 1+1-dimensional d'Alembertian:

(define (make-wave-operator c)
  (lambda (nodes)
    (make-operator
     (lambda (f)
       (vector (d/dt f) (basis:* (- c) (d/dx f))))
     (lambda (f)
       (vector (d/dt f) (basis:* c (d/dx f))))
     (lambda (v w)
       (basis:* -1 (basis:dot v w))))))


;;; Characteristic bending: For t < t1, the operator agrees with the wave
;;; operator.  For t > t2, the equation becomes elliptic.  The characteristics
;;; are "bent" between t1 and t2.

(define (make-bent-operator c t1 t2)
  (let ((phi (make-bending-coeff t1 t2)))
    (lambda (nodes)
      (let ((phi (function->poly phi nodes)))
	(make-operator
	 (lambda (f)
	   (vector (basis:* phi (d/dt f)) (basis:* (- (square c)) (d/dx f))))
	 (lambda (f)
	   (vector (d/dt f) (d/dx f)))
	 (lambda (v w)
	   (basis:* -1 (basis:dot v w))))))))


;;; Let this be a polynomial for now:

(define (cut-off t)

  ;; This one lets the wave operator transition nicely into a "parabolic"
  ;; operator.  (But without a time derivative!)

  (cond ((<= t 0) 1.)
	((>= t 1) 0.)
	(else (+ (* 2 (cube t)) (* -3 (square t)) 1))))

(define (cut-off-1 t)

  ;; This one lets the wave operator transition into an elliptic operator.

  (if (> t 0)
      (- 1 (cube t))
      1.))

(define (make-bending-coeff t1 t2)
  (let ((delta (- t2 t1)))
    (lambda (node)
      (let ((t (node:get-y node)))
	(cut-off (/ (- t t1) delta))))))


;;; The Laplcian operator for "real" functions:

(define (real-laplacian nodes)
  (make-operator
   real-gradient
   real-gradient
   (lambda (v w) (basis:scalar* -1 (basis:dot v w)))))
