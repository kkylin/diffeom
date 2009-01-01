;;; This file uses ScmUtils to derive the Lagrangian for rigid bodies using an
;;; Euler-like coordinate system that covers a region of SO(3) different from
;;; standard Euler angles.


;;; These are the Euler-like charts we are actually working with:


;;; A canonical rotation matrix that isn't defined in ScmUtils:

(define (rotate-y angle)
  (vector
   (vector (cos angle)     0           (- (sin angle)))
   (vector      0          1                   0)
   (vector (sin angle)     0           (cos angle))))

;;; Some procedures for testing charts built on Euler angles:

(define (compare-euler-angles chart euler->m)
  (let ((theta (random pi))
	(phi (- (random (* 2 pi)) pi))
	(psi (- (random (* 2 pi)) pi)))
    (let ((A (chart:coords->point (vector psi theta phi) chart))
	  (B (euler->m (vector theta phi psi))))
    (print-matrix A)
    (print-matrix B)
    (write-line `(diff = ,(matrix:max (matrix:- A B)))))))

(define (test-euler-chart chart)
  (let ((psi (- (random (* 2 pi)) pi))
	(theta (random pi))
	(phi (- (random (* 2 pi)) pi)))
    (let ((v (vector psi theta phi)))
      (vector:distance (chart:point->coords (chart:coords->point v chart)
					    chart)
		       v))))

(define (test-euler-tangent-chart chart range)
  (let ((psi (- (random (* 2 pi)) pi))
	(theta (random pi))
	(phi (- (random (* 2 pi)) pi))
	(psidot (- (random (* 2 range)) range))
	(thetadot (- (random (* 2 range)) range))
	(phidot (- (random (* 2 range)) range))

	(chart (make-tangent-chart chart)))
    (let ((v (vector psi theta phi psidot thetadot phidot)))
      (vector:distance (chart:point->coords (chart:coords->point v chart)
					    chart)
		       v))))

(define (test-euler-tt-chart chart range)
  (let ((psi (- (random (* 2 pi)) pi))
	(theta (random pi))
	(phi (- (random (* 2 pi)) pi))
	(psidot (- (random (* 2 range)) range))
	(thetadot (- (random (* 2 range)) range))
	(phidot (- (random (* 2 range)) range))
	(a (- (random (* 2 range)) range))
	(b (- (random (* 2 range)) range))
	(c (- (random (* 2 range)) range))
	(d (- (random (* 2 range)) range))
	(e (- (random (* 2 range)) range))
	(f (- (random (* 2 range)) range))

	(chart (make-tangent-chart (make-tangent-chart chart))))
    (let ((v (vector psi theta phi psidot thetadot phidot
		     a b c d e f)))
      (vector:distance (chart:point->coords (chart:coords->point v chart)
					    chart)
		       v))))

;;; Generate rotation matrix from angles: (Why is the order of arguments used
;;; in the book so *odd*?)

(define (euler0->m angles)
  (let ((theta (vector-ref angles 0))
	(phi (vector-ref angles 1))
	(psi (vector-ref angles 2)))
    (matrix:* (rotate-z phi)
	      (rotate-y (- theta))
	      (rotate-z psi))))

(define euler1->m
  (let ((rot (vector (vector -1 0 0)
		     (vector 0 -1 0)
		     (vector 0 0 1))))
    (lambda (angles)
      (let ((theta (vector-ref angles 0))
	    (phi (vector-ref angles 1))
	    (psi (vector-ref angles 2)))
	(matrix:* (rotate-z phi)
		  (rotate-y (- theta))
		  (rotate-z psi)
		  rot)))))

(define euler2->m
  (let ((rot-x (vector (vector 1  0  0)
		       (vector 0  0  1)
		       (vector 0 -1  0)))
	(rot-y (vector (vector -1  0  0)
		       (vector  0  1  0)
		       (vector  0  0 -1))))
    (lambda (angles)
      (let ((theta (vector-ref angles 0))
	    (phi (vector-ref angles 1))
	    (psi (vector-ref angles 2)))
	(matrix:* (rotate-y phi)
		  rot-y
		  (rotate-z (- theta))
		  rot-x
		  (rotate-z psi))))))

(define euler3->m
  (let ((rot-x (vector (vector 1  0  0)
		       (vector 0  0  1)
		       (vector 0 -1  0)))
	(rot-y (vector (vector -1  0  0)
		       (vector  0  1  0)
		       (vector  0  0 -1)))
	(rot-z (vector (vector -1  0  0)
		       (vector  0 -1  0)
		       (vector  0  0  1))))
    (lambda (angles)
      (let ((theta (vector-ref angles 0))
	    (phi (vector-ref angles 1))
	    (psi (vector-ref angles 2)))
	(matrix:* (rotate-y phi)
		  rot-y
		  (rotate-z (- theta))
		  rot-x
		  (rotate-z psi)
		  rot-z)))))


;;; Generate the angular velocity vector:

(define (((make-euler->omega euler->m) angles-path) t)
  (define (m-on-path t)
    (euler->m (angles-path t)))
  (define (w-cross t)
    (matrix:* ((derivative m-on-path) t)
	      (matrix:transpose (m-on-path t))))
  (antisymmetric->3vector-components (w-cross t)))


;;; Generate the angular velocity vector in the body frame:

(define (((make-euler->omega-body euler->m) angles-path) t)
  (matrix:matrix*vector
   (matrix:transpose (euler->m (angles-path t)))
   (((make-euler->omega euler->m) angles-path) t)))


;;; Compute the expression for the angular velocity in terms of the angles:

(define (angular-velocity-expression euler->m)
  (let ((euler->omega-body (make-euler->omega-body euler->m)))
    (show-time
     (lambda ()
       ((compose ham:simplify easy-simplify)
	((euler->omega-body
	  (vector (literal-function 'theta)
		  (literal-function 'phi)
		  (literal-function 'psi)))
	 't))))))


;;; Compute the kinetic energy!

(define ((t-rigid-body-0 a b c) state)
  (let ((q (state->q state))
	(qdot (state->qdot state))
	(t (state->t state)))
    (let ((theta (vector-ref q 0))
	  (phi (vector-ref q 1))
	  (psi (vector-ref q 2))
	  (thetadot (vector-ref qdot 0))
	  (phidot (vector-ref qdot 1))
	  (psidot (vector-ref qdot 2)))
      (let ((w_a (- (* thetadot (sin psi)) (* (sin theta) (cos psi) phidot)))
	    (w_b (+ (* (cos psi) thetadot) (* (sin theta) (sin psi) phidot)))
	    (w_c (+ psidot (* (cos theta) phidot))))
	(* 1/2 (+ (* a (square w_a))
		  (* b (square w_b))
		  (* c (square w_c))))))))

(define ((t-rigid-body-1 a b c) state)
  (let ((q (state->q state))
	(qdot (state->qdot state))
	(t (state->t state)))
    (let ((theta (vector-ref q 0))
	  (phi (vector-ref q 1))
	  (psi (vector-ref q 2))
	  (thetadot (vector-ref qdot 0))
	  (phidot (vector-ref qdot 1))
	  (psidot (vector-ref qdot 2)))
      (let ((w_a (- (* (sin theta) (cos psi) phidot) (* thetadot (sin psi))))
	    (w_b (- (+ (* (sin theta) (sin psi) phidot)
		       (* (cos psi) thetadot))))
	    (w_c (+ psidot (* (cos theta) phidot))))
	(* 1/2 (+ (* a (square w_a))
		  (* b (square w_b))
		  (* c (square w_c))))))))

(define ((t-rigid-body-2 a b c) state)
  (let ((q (state->q state))
	(qdot (state->qdot state))
	(t (state->t state)))
    (let ((theta (vector-ref q 0))
	  (phi (vector-ref q 1))
	  (psi (vector-ref q 2))
	  (thetadot (vector-ref qdot 0))
	  (phidot (vector-ref qdot 1))
	  (psidot (vector-ref qdot 2)))
      (let ((w_a (+ (* (sin theta) (cos psi) phidot) (* (sin psi) thetadot)))
	    (w_b (- (* (cos psi) thetadot) (* (sin theta) (sin psi) phidot)))
	    (w_c (- psidot (* (cos theta) phidot))))
	(* 1/2 (+ (* a (square w_a))
		  (* b (square w_b))
		  (* c (square w_c))))))))

(define ((t-rigid-body-3 a b c) state)
  (let ((q (state->q state))
	(qdot (state->qdot state))
	(t (state->t state)))
    (let ((theta (vector-ref q 0))
	  (phi (vector-ref q 1))
	  (psi (vector-ref q 2))
	  (thetadot (vector-ref qdot 0))
	  (phidot (vector-ref qdot 1))
	  (psidot (vector-ref qdot 2)))
      (let ((w_a (- (+ (* (sin theta) (cos psi) phidot)
		       (* (sin psi) thetadot))))
	    (w_b (- (* (sin theta) (sin psi) phidot) (* (cos psi) thetadot)))
	    (w_c (- psidot (* (cos theta) phidot))))
	(* 1/2 (+ (* a (square w_a))
		  (* b (square w_b))
		  (* c (square w_c))))))))


;;; Compute the angular momentum!

(define ((state->L-body-0 a b c) state)
  (let ((q (state->q state))
	(qdot (state->qdot state))
	(t (state->t state)))
    (let ((theta (vector-ref q 0))
	  (phi (vector-ref q 1))
	  (psi (vector-ref q 2))
	  (thetadot (vector-ref qdot 0))
	  (phidot (vector-ref qdot 1))
	  (psidot (vector-ref qdot 2)))
      (let ((w_a (- (* thetadot (sin psi)) (* (sin theta) (cos psi) phidot)))
	    (w_b (+ (* (cos psi) thetadot) (* (sin theta) (sin psi) phidot)))
	    (w_c (+ psidot (* (cos theta) phidot))))
	(vector (* a w_a) (* b w_b) (* c w_c))))))

(define ((state->L-body-1 a b c) state)
  (let ((q (state->q state))
	(qdot (state->qdot state))
	(t (state->t state)))
    (let ((theta (vector-ref q 0))
	  (phi (vector-ref q 1))
	  (psi (vector-ref q 2))
	  (thetadot (vector-ref qdot 0))
	  (phidot (vector-ref qdot 1))
	  (psidot (vector-ref qdot 2)))
      (let ((w_a (- (* (sin theta) (cos psi) phidot) (* thetadot (sin psi))))
	    (w_b (- (+ (* (sin theta) (sin psi) phidot)
		       (* (cos psi) thetadot))))
	    (w_c (+ psidot (* (cos theta) phidot))))
	(vector (* a w_a) (* b w_b) (* c w_c))))))

(define ((state->L-body-2 a b c) state)
  (let ((q (state->q state))
	(qdot (state->qdot state))
	(t (state->t state)))
    (let ((theta (vector-ref q 0))
	  (phi (vector-ref q 1))
	  (psi (vector-ref q 2))
	  (thetadot (vector-ref qdot 0))
	  (phidot (vector-ref qdot 1))
	  (psidot (vector-ref qdot 2)))
      (let ((w_a (+ (* (sin theta) (cos psi) phidot) (* (sin psi) thetadot)))
	    (w_b (- (* (cos psi) thetadot) (* (sin theta) (sin psi) phidot)))
	    (w_c (- psidot (* (cos theta) phidot))))
	(vector (* a w_a) (* b w_b) (* c w_c))))))

(define ((state->L-body-3 a b c) state)
  (let ((q (state->q state))
	(qdot (state->qdot state))
	(t (state->t state)))
    (let ((theta (vector-ref q 0))
	  (phi (vector-ref q 1))
	  (psi (vector-ref q 2))
	  (thetadot (vector-ref qdot 0))
	  (phidot (vector-ref qdot 1))
	  (psidot (vector-ref qdot 2)))
      (let ((w_a (- (+ (* (sin theta) (cos psi) phidot)
		       (* (sin psi) thetadot))))
	    (w_b (- (* (sin theta) (sin psi) phidot) (* (cos psi) thetadot)))
	    (w_c (- psidot (* (cos theta) phidot))))
	(vector (* a w_a) (* b w_b) (* c w_c))))))

(define ((state->L-space-0 a b c) state)
  (let ((angles (state->q state)))
    (* (euler0->m angles) ((state->L-body-0 a b c) state))))

(define ((state->L-space-1 a b c) state)
  (let ((angles (state->q state)))
    (* (euler1->m angles) ((state->L-body-1 a b c) state))))

(define ((state->L-space-2 a b c) state)
  (let ((angles (state->q state)))
    (* (euler2->m angles) ((state->L-body-2 a b c) state))))

(define ((state->L-space-3 a b c) state)
  (let ((angles (state->q state)))
    (* (euler3->m angles) ((state->L-body-3 a b c) state))))


;;; The state derivatives:

(define (rigid-sysder-0 a b c)
  (lagrangian->state-derivative
   (t-rigid-body-0 a b c)))

(define (rigid-sysder-1 a b c)
  (lagrangian->state-derivative
   (t-rigid-body-1 a b c)))

(define (rigid-sysder-2 a b c)
  (lagrangian->state-derivative
   (t-rigid-body-2 a b c)))

(define (rigid-sysder-3 a b c)
  (lagrangian->state-derivative
   (t-rigid-body-3 a b c)))


;;; This is rather useful:

(define rigid-qqdot
  (->state 't
	   (vector 'theta 'phi 'psi)
	   (vector 'thetadot 'phidot 'psidot)))

(define rigid-qp
  (->state 't
	   (vector 'theta 'phi 'psi)
	   (vector 'p_theta 'p_phi 'p_psi)))


;;; Something useful to remember (thanks to CPH): This should fix the
;;; simplifier so that the really big expressions can be simplified.

;(ge '(user))
;(in-package scmutils-base-environment
;  ((pcf-package 'set-gcd-method!) (pcf-package 'gcd-euclid)))
;(ge generic-environment)


;;; Stupid exchange of order of arguments:

(define (traditional->correct-order v)
  (let ((theta (vector-ref v 0))
	(phi (vector-ref v 1))
	(psi (vector-ref v 2))

	(thetadot (vector-ref v 3))
	(phidot (vector-ref v 4))
	(psidot (vector-ref v 5)))
    (vector psi theta phi psidot thetadot phidot)))
