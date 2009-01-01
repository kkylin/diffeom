(load "load-ode")

(define (make-v.field a b c)
  (lambda (state)
    (vector-append (vector 1)
		   (rigid-field-1 a b c
				  (state->q state)
				  (state->qdot state)))))

(define (reorder v)
  (let ((t (state->t v))
	(q (state->q v))
	(qdot (state->qdot v)))
    (let ((psi (vector-ref q 0))
	  (theta (vector-ref q 1))
	  (phi (vector-ref q 2))
	  (psidot (vector-ref qdot 0))
	  (thetadot (vector-ref qdot 1))
	  (phidot (vector-ref qdot 2)))
      (->state t (vector theta phi psi) (vector thetadot phidot psidot)))))

(define v.field (make-v.field 1. (sqrt 2) 2.))
(define step (compose integrator:get-new-x (make-rk4-integrator .01)))

(define results
  (show-time
   (lambda ()
     (let loop ((i 0)
		(x (->state 0 (vector 0 1 0) (vector .1 .1 .1)))
		(results '()))
       (if (< i 10000)
	   (loop (+ i 1) (step x v.field (lambda () 'foo)) (cons x results))
	   (map reorder
		(sort results
		      (lambda (x y) (< (state->t x) (state->t y))))))))))

(define E-errors (map (t-rigid-body 1. (sqrt 2) 2.) results))
(define L-errors (map (state->L-space 1. (sqrt 2) 2.) results))
(define L1-errors (map vector-first L-errors))
(define L2-errors (map vector-second L-errors))
(define L3-errors (map vector-third L-errors))

;;; Graphics devices:

(define dev 'undefined)

(define (open)
  (if (eq? dev 'undefined)
      (begin
	(set! dev (make-graphics-device 'x))
	(graphics-operation dev 'set-background-color "white")
	(graphics-operation dev 'set-foreground-color "blue")
	(graphics-operation dev 'set-mouse-color "black")
	(graphics-set-coordinate-limits dev -2 -2 2 2)
	(graphics-clear dev))))

(define (close)
  (if (not (eq? dev 'undefined))
      (begin
	(graphics-close dev)
	(set! dev 'undefined))))

(define (clear)
  (open)
  (graphics-clear dev))

(define (plot-conservation-error l)
  (let* ((ref (car l))
	 (l (map (lambda (x) (relative-error x ref)) l))
	 (max (apply max l))
	 (min (apply min l))
	 (len (- (length l) 1))
	 (ref (car l)))
    (write-line `(range: ,min to ,max))
    (open)
    (graphics-enable-buffering dev)
    (graphics-set-coordinate-limits dev 0. min len max)
    (graphics-move-cursor dev 0. 0.)
    (let loop ((i 1) (l (cdr l)))
      (if (null? l)
	  (graphics-disable-buffering dev)
	  (begin
	    (graphics-drag-cursor dev i (car l))
	    (loop (+ i 1) (cdr l)))))))

(define (set-fg color)
  (graphics-operation dev 'set-foreground-color color))
