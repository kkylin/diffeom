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

;;; This file defines basis functions that are still polynomial, but are
;;; represented by real Scheme procedures and can thus undergo general
;;; coordinate transformations in a nice way.  This is not so important here
;;; (in fact, it is a slower and less accurate implementation), but is useful
;;; for extending FEM to manifolds.

(declare (usual-integrations))


;;; Constructor:

(define (proc->real f)
  (package-basis-function-methods
   '2d-real-basis-function
   f
   f
   (make-real-adder f)
   (make-real-subtractor f)
   (make-real-multiplier f)
   (make-real-scalar-multiplier f)))

(define make-real-basis-function
  (compose proc->real basis-function->function make-polynomial-basis-function))


;;; Operations on basis functions:

(define (make-real-adder f)
  (lambda (g)
    (let ((g (basis:get-rep g)))
      (proc->real
       (lambda (x)
	 (+ (f x) (g x)))))))

(define (make-real-subtractor f)
  (lambda (g) 
    (let ((g (basis:get-rep g)))
      (proc->real
       (lambda (x)
	 (- (f x) (g x)))))))

(define (make-real-multiplier f)
  (lambda (g)
    (let ((g (basis:get-rep g)))
      (proc->real
       (lambda (x)
	 (* (f x) (g x)))))))

(define (make-real-scalar-multiplier f)
  (lambda (a)
    (proc->real
     (lambda (x)
       (* a (f x))))))


;;; The gradient is needed for defining the laplacian:

(define (real-gradient f)
  (let* ((f (compose vector (basis:get-rep f)))
	 (fx (proc->real (compose vector-first ((pdiff 0) f))))
	 (fy (proc->real (compose vector-first ((pdiff 1) f)))))
    (vector fx fy)))
