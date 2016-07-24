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

;;; This file defines some structures related to vector fields:


;;; Here's a trivial vector field on the circle:

(define (circle-field p)
  (let ((x (vector-ref p 0))
	(y (vector-ref p 1)))
    (imbedding->tangent circle p (vector (- y) x))))


;;; Here is the nonlinear pendulum.  This shouldn't be *that* hard to define,
;;; but it is.  Why?  What should we change about the system?

;;; Related to the problem of defining vector fields is the issue of
;;; efficiency.  The definition problem can be solved by making more charts,
;;; but efficiency would suffer even more.  How do we fix that?

(define (make-pendulum g mass length)
  (let* ((k (/ g length))
	 (-k (* -1 k))
	 (chart-1 (make-spherical-chart 1 '(0 1) 0))
	 (chart-2 (make-spherical-chart 1 '(1 0) 0))
	 (find-chart
	  (make-tangent-chart-finder (lambda (x)
				       (if (chart:member? x chart-1)
					   chart-1
					   chart-2)))))

    (lambda (p)

      ;; p should be a point from the cylinder constructed above.

      (let* ((x (tangent:get-anchor p))
	     (chart (if (chart:member? x chart-1) chart-1 chart-2)))

	(let ((xdot (chart:push-forward p chart)))
	  (make-tangent (find-chart p)
			p
			(vector-append
			 xdot
			 (vector (* (if (chart:member? x chart-1) -k k)
				    (vector-ref x 1))))))))))


;;; This is useful for checking how well the integrator is doing:

(define (make-pendulum-energy-function g mass length)
  (lambda (p)
    (let ((x (tangent:get-anchor p))
	  (v (tangent->imbedded-velocity circle p)))
      (- (/ (* mass (vector:magnitude^2 v)) 2)
	 (* mass g length (vector-ref x 0))))))


;;; Spherical pendulum:

(define make-spherical-pendulum
  (let* ((C1 (make-tangent-chart (make-spherical-chart 2 '(2 0 1) 0)))
	 (C2 (make-tangent-chart (make-spherical-chart 2 '(1 0 2) pi)))
	 (TS^2 (charts->manifold (list C1 C2))))
    (lambda (g mass length)
      (let* ((k (/ g length))
	     (-k (- k)))
	(lambda (p)
	  (let* ((chart (manifold:find-best-chart TS^2 p))
		 (x (chart:point->coords p chart))
		 (phi (vector-ref x 0))
		 (theta (vector-ref x 1))
		 (phidot (vector-ref x 2))
		 (thetadot (vector-ref x 3)))
	    (make-tangent chart
			  p

			  ;; Ended up deriving these things using Lagrangian
			  ;; mechanics anyway; might as well automate it.

			  (if (eq? chart C1)
			      (vector phidot
				      thetadot
				      (* (+ 1 (* (cos phi) (square thetadot)))
					 (sin phi))
				      (* -2 (cot phi) phidot thetadot))
			      (vector phidot
				      thetadot
				      (* (+ (* (sin phi) (square thetadot))
					    (sin theta))
					 (cos phi))
				      (+ (* -2 (cot phi) phidot thetadot)
					 (/ (cos theta) (sin phi))))))))))))

;;; And in phase space:

(define make-spherical-H-pendulum
  (let* ((C1 (make-cotangent-chart (make-spherical-chart 2 '(2 0 1) 0)))
	 (C2 (make-cotangent-chart (make-spherical-chart 2 '(1 0 2) pi)))
	 (T*S^2 (charts->manifold (list C1 C2))))
    (lambda (g mass length)
      (let ((k1 (/ (* mass (square length))))
	    (k2 (* mass g length)))
	(lambda (p)
	  (let* ((chart (manifold:find-best-chart T*S^2 p))
		 (x (chart:point->coords p chart))
		 (phi (vector-ref x 0))
		 (theta (vector-ref x 1))
		 (p_phi (vector-ref x 2))
		 (p_theta (vector-ref x 3)))
	    (make-tangent chart
			  p
			  (if (eq? chart C1)
			      (vector (* k1 p_phi)
				      (* (/ k1 (square (sin phi))) p_theta)
				      (+ (* k1 (square p_theta)
					    (/ (* (square (sin phi))
						  (tan phi))))
					 (* k2 (sin phi)))
				      0)
			      (vector (* k1 p_phi)
				      (* (/ k1 (square (sin phi))) p_theta)
				      (+ (* k1 (square p_theta)
					    (/ (* (square (sin phi))
						  (tan phi))))
					 (* k2 (cos phi) (sin theta)))
				      (* k2 (sin phi) (cos theta)))))))))))


;;; An example of something defined using Lagrangian methods:

(define R^3 (make-euclidean-space 3))
(define TR^3 (make-tangent-bundle R^3))

(define (make-free-particle-lagrangian m)
  (let ((m/2 (/ m 2)))
    (make-real-map TR^3
		   (lambda (p)
		     (* m/2 (vector:magnitude^2 (tangent:get-coords p)))))))


;;; An example of something defined using Hamiltonian methods:

(define T*R^3 (make-cotangent-bundle R^3))

(define (make-free-particle-hamiltonian m)
  (let ((m/2 (/ m 2)))
    (make-real-map T*R^3
		   (lambda (p)
		     (* m/2 (vector:magnitude^2 (cotangent:get-coords p)))))))

;;; The Lagrangian for a rather simple potential:

(define TR^3 (make-tangent-bundle R^3))

(define falling-lagrangian
  (make-real-map
   TR^3 (lambda (p)
	  (- (* 1/2 (vector:magnitude^2 (tangent:get-coords p)))
	     (vector-third (tangent:get-anchor p))))))


;;; And the equivalent Hamiltonian:

(define T*R^3 (make-cotangent-bundle R^3))

(define falling-hamiltonian
  (make-real-map
   T*R^3 (lambda (p)
	   (+ (* 1/2 (vector:magnitude^2 (cotangent:get-coords p)))
	      (vector-third (cotangent:get-anchor p))))))


;;; And now for rigid bodies:

(define SO3 (make-rotational-group))
(define TSO3 (make-tangent-bundle SO3))
(define T*SO3 (make-cotangent-bundle SO3))

;;; What is the easiest way to make a Lagrange top?  Just write down the
;;; Lagrangian!

(define (antisymmetric-matrix->vector A)
  (vector (matrix-ref A 2 1) (matrix-ref A 0 2) (matrix-ref A 1 0)))

(define (tangent->angular-velocity p)
  (let* ((M (tangent:get-anchor p))
	 (chart (tangent:get-chart p))
	 (Mdot
	  (vector->matrix 3 3
			  (push-forward-in-coords
			   (compose matrix:flatten
				    (chart:get-inverse-map chart))
			   (chart:point->coords M chart)
			   (tangent:get-coords p)))))
    (antisymmetric-matrix->vector (matrix:* Mdot (transpose M)))))

(define (make-free-rigid-body-angular-momentum A B C)
  (make-simple-map
   TSO3
   (make-euclidean-space 3)
   (lambda (p)
     (let* ((M (tangent:get-anchor p))
	    (w (apply-linear-transformation
		(transpose M) (tangent->angular-velocity p)))
	    (w-prime (apply-linear-transformation (transpose M) w)))
       (vector (* A (vector-ref w-prime 0))
	       (* B (vector-ref w-prime 1))
	       (* C (vector-ref w-prime 2)))))))

(define (make-free-rigid-body-lagrangian A B C)
  (make-real-map
   TSO3
   (lambda (p)
     (let* ((M (tangent:get-anchor p))
	    (w (apply-linear-transformation
		(transpose M) (tangent->angular-velocity p)))
	    (w-prime (apply-linear-transformation (transpose M) w)))
       (* 1/2
	  (+ (* A (square (vector-ref w-prime 0)))
	     (* B (square (vector-ref w-prime 1)))
	     (* C (square (vector-ref w-prime 2)))))))))

(define (free-rigid-body-field-maker a b c)
  (let* ((charts (manifold:get-finite-atlas TSO3))
	 (fields (list rigid-field-0 rigid-field-1
		       rigid-field-2 rigid-field-3))
	 (charts&fields (map cons charts fields)))
    (lambda (chart)
      (let ((field (assq chart charts&fields)))
	(if field
	    (lambda (x)
	      ((cdr field) a b c (vector-head x 3) (vector-tail x 3)))
	    (error "Unknown chart! -- FREE-RIGID-BODY"))))))

(define (scmutils-rigid-body-field-maker a b c)
  (let* ((charts (manifold:get-finite-atlas TSO3))
	 (fields (map (make-sysder a b c)
		      (list t-rigid-body-0
			    t-rigid-body-1
			    t-rigid-body-2
			    t-rigid-body-3)))
	 (charts&fields (map cons charts fields)))
    (lambda (chart)
      (let ((field (assq chart charts&fields)))
	(if field
	    (cdr field)
	    (error "Unknown chart! -- FREE-RIGID-BODY"))))))

(define (make-sysder a b c)
  (lambda (t-rigid-body)
    (show-time
     (lambda ()
       (let* ((sysder (lambda (a b c)
			(lagrangian->state-derivative
			 (t-rigid-body a b c))))
	      (compiled (compile-sysder 3 sysder))
	      (field (compiled a b c)))
	 (lambda (x)
	   (state->manifold (field (manifold->state x)))))))))

(define (make-rigid-body-energy a b c)
  (let* ((charts (cons euler-angles (manifold:get-finite-atlas SO3)))
	 (energies (list (t-rigid-body a b c)
			 (t-rigid-body-0 a b c)
			 (t-rigid-body-1 a b c)
			 (t-rigid-body-2 a b c)
			 (t-rigid-body-3 a b c)))
	 (charts&energies (map cons charts energies)))
    (lambda (tangent)
      (let* ((chart (tangent:get-chart tangent))
	     (energy (assq chart charts&energies)))
	(if energy
	    ((cdr energy)
	     (manifold->state
	      (vector-append
	       (chart:point->coords (tangent:get-anchor tangent) chart)
	       (tangent:get-coords tangent))))
	    (error "Unknown chart! -- FREE-RIGID-BODY"))))))

(define (make-rigid-body-momentum a b c)

  ;; This version is more accurate than MAKE-RIGID-BODY-ANGULAR-MOMENTUM.

  (let* ((charts (cons euler-angles (manifold:get-finite-atlas SO3)))
	 (momenta (list (state->L-space a b c)
			 (state->L-space-0 a b c)
			 (state->L-space-1 a b c)
			 (state->L-space-2 a b c)
			 (state->L-space-3 a b c)))
	 (charts&momenta (map cons charts momenta)))
    (lambda (tangent)
      (let* ((chart (tangent:get-chart tangent))
	     (momentum (assq chart charts&momenta)))
	(if momentum
	    ((cdr momentum)
	     (manifold->state
	      (vector-append
	       (chart:point->coords (tangent:get-anchor tangent) chart)
	       (tangent:get-coords tangent))))
	    (error "Unknown chart! -- FREE-RIGID-BODY"))))))

(define (manifold->state x)
  (let ((psi (vector-ref x 0))
	(theta (vector-ref x 1))
	(phi (vector-ref x 2))
	(psidot (vector-ref x 3))
	(thetadot (vector-ref x 4))
	(phidot (vector-ref x 5)))
    (->state 0 (vector theta phi psi) (vector thetadot phidot psidot))))

(define (state->manifold state)
  (let ((q (state->q state))
	(qdot (state->qdot state))
	(t (state->t state)))
    (let ((theta (vector-ref q 0))
	  (phi (vector-ref q 1))
	  (psi (vector-ref q 2))

	  (thetadot (vector-ref qdot 0))
	  (phidot (vector-ref qdot 1))
	  (psidot (vector-ref qdot 2)))
      (vector psi theta phi psidot thetadot phidot))))
