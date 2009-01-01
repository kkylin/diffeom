;;; This file uses ScmUtils to derive the Lagrangian for tops.  This is much
;;; like the rigid body stuff.

(load "rigid")


;;; The Lagrangian for the top (note that the moments of inertia are with
;;; respect to the pivot):

(define (make-top-lagrangian-1 a b c MgR)
  (let ((T (T-rigid-body a b c)))
    (lambda (state)
      (let ((theta (vector-ref (state->q state) 0)))
	(- (T state) (* MgR (cos theta)))))))

(define (make-top-lagrangian-2 a b c MgR)
  (let ((T (kT-rigid-body a b c)))
    (lambda (state)
      (let* ((q (state->q state))
	     (theta (vector-ref q 0))
	     (phi (vector-ref q 1)))
	(+ (T state) (* MgR (cos phi) (sin theta)))))))
	

;;; The corresponding Euler-Lagrange equations:

(define (top-sysder-1 a b c MgR)
  (lagrangian->state-derivative
   (make-top-lagrangian-1 a b c MgR)))

(define (top-sysder-2 a b c MgR)
  (lagrangian->state-derivative
   (make-top-lagrangian-2 a b c MgR)))


;;; Compute some expressions:

(let ((port (open-output-file "foo")))
  (pp (traditional->correct-order
       (vector-tail
	(show-time
	 (lambda ()
	   (*sysder-simplify*
	    ((top-sysder-1 'a 'b 'c 'MgR) rigid-qqdot))))
	1))
      port)
  (close-output-port port))

(let ((port (open-output-file "bar")))
  (pp (traditional->correct-order
       (vector-tail
	(show-time
	 (lambda ()
	   (*sysder-simplify*
	    ((top-sysder-2 'a 'b 'c 'MgR) rigid-qqdot))))
	1))
      port)
  (close-output-port port))
