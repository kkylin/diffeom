(load "load-ode")

(define result
  (show-time
   (lambda ()
     (rigid-body-path singular-init 1.))))

(define e-list
  (show-time
   (lambda ()
     (map (compose vector-first (make-rigid-body-energy 1. (sqrt 2) 2.) cadr)
	  result))))

(define e-errors
  (let ((ref (car e-list)))
    (show-time
     (lambda ()
       (map (lambda (val) (relative-error val ref)) e-list)))))
