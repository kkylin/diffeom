;;; This uses ScmUtils to compute the usual axisymmetric top:

(set! *ode-integration-method*  'bulirsch-stoer)

(define (my-ode-advancer v x dt tol)
  (let ((dt/2 (/ dt 2.))
	(dt/6 (/ dt 6.)))
    (let* ((F1 (v x))
	   (F2 (v (vector:+ x (vector:scalar*vector dt/2 F1))))
	   (F3 (v (vector:+ x (vector:scalar*vector dt/2 F2))))
	   (F4 (v (vector:+ x (vector:scalar*vector dt F3)))))
      (vector:+ (vector:scalar*vector dt/6
				      (vector:+ F1
						(vector:scalar*vector 2. F2)
						(vector:scalar*vector 2. F3)
						F4))
		x))))

(define (rigid-body-evolver a b c x0 t-final dt tol)
  (let ((v.field (compiled-rigid-body a b c)))
    (let loop ((x x0) (results '()))
      (if (< (state->t x) t-final)
	  (let ((new-x (ode-advancer v.field x dt tol)))
	    (loop new-x (cons x results)))
	  results))))


;;; The vector field:

(define (rigid-body-sysder a b c)
  (lagrangian->state-derivative
   (t-rigid-body a b c)))

(newline)
(display "*** Compiling state derivative")

(define compiled-rigid-body
  (show-time
   (lambda ()
     (compile-sysder 3 rigid-body-sysder))))


;;; Useful functions:

;((STATE->L-SPACE A B C) STATE) => angular momentum in the reference frame.
;((STATE->L-BODY A B C) STATE) => angular momentum in the body frame.
;(RELATIVE-ERROR VALUE REFERENCE-VALUE) => error of VALUE relative to
;   REFERENCE-VALUE.


;;; Initial conditions from the book:

;;; In Euler coordinates:
;;; q0 = #(1 0 0), qdot0 = #(0.1 0.1 0.1).
;;; Step size = 0.01, and final time is 100.0.
;;; Maximum local truncation error is 1.0e-12.
;;; A=1, B=sqrt(2), C=2.

;;; Note that the order of components needs to be switched when using these
;;; initial conditions with the manifold stuff.

(newline)
(display "*** Evolving trajectories")

(define result
  (show-time
   (lambda ()
     (rigid-body-evolver 1 (sqrt 2) 2
			 (->state 0. (vector 1. 0. 0.) (vector -.1 -.01 -.01))
			 100.
			 .01 1.0e-12))))

(newline)
(display "*** Saving results")

(let ((port (open-output-file "rigid-reg.data")))
  (for-each
   (lambda (state)
     (display (state->t state) port)
     (display " " port)
     (display (state->q state) port)
     (display " " port)
     (display (state->qdot state) port)
     (newline port))
   (sort results (lambda (x y) (< (state->t x) (state->t y)))))
  (close-output-port port))

;;; Directly from the text:

(define (do-it A B C state0 final-t dt tol)
  (let ((dstate (compiled-rigid-body A B C))
	(L0 ((state->L-space A B C) state0))
	(E0 ((T-rigid-body A B C) state0)))
    (let ((Lx0 (vector-ref L0 0))
	  (Ly0 (vector-ref L0 1))
	  (Lz0 (vector-ref L0 2)))
      (let loop ((state state0))
	(if (< (state->t state) final-t)
	    (let ((ns (ode-advancer dstate state dt tol)))
	      (let ((L ((state->L-space A B C) ns))
		    (E ((T-rigid-body A B C) ns))
		    (t (state->t ns)))
		(let ((Lx (vector-ref L 0))
		      (Ly (vector-ref L 1))
		      (Lz (vector-ref L 2)))
		  (let ((error-Lx (relative-error Lx Lx0))
			(error-Ly (relative-error Ly Ly0))
			(error-Lz (relative-error Lz Lz0))
			(error-E (relative-error E E0)))
		    (plot-point window t error-Lx)
		    (plot-point window t error-Ly)
		    (plot-point window t error-Lz)
		    (plot-point window t error-E)
		    (loop ns))))))))))

#|
;;; Comments by GJS:

;;; For QC Runge-Kutta 4 
(set! *ode-integration-method*  (quote qcrk4))
(define window (frame 0. 100. -1.e-12 1.e-12))

;;; For bulirsch-stoer 
(set! *ode-integration-method*  (quote bulirsch-stoer))
(define window (frame 0. 100. -1.e-13 1.e-13))

;;; Comes by the coordinate singularity several times

(do-it 1. (sqrt 2.) 2.
       (->state 0.0
		(vector 1. 0. 0.)
		(vector 0.1 0.1 0.1))
       100.0
       0.01
       1.0e-12)


;;; Whizzing rather close to a singularity

(do-it 1. (sqrt 2.) 2.
       (->state 0.0
		(vector 1. 0. 0.)
		(vector -0.1 -.01 -.01))
       100.0
       0.01
       1.0e-12)
|#
