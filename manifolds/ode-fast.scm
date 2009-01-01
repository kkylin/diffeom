;;; This is like ode.scm, only the integrator has its distortion checks
;;; disabled.

(declare (usual-integrations))

(define (fast-v.field->flow M make-local-field next-step)
  (lambda (p t-final . aux)

    ;; Reset the error-reporting mechanism:
    (set! *point-of-failure* #f)

    ;; AUX lets the user specify an initial time (optional).

    (let next-point ((p p)
		     (t (if (not (null? aux)) (car aux) 0.))
		     (result '()))
      (if (< t t-final)

	  ;; Find all the charts containing this point and try each of them:

	  (let ((charts (manifold:get-local-atlas M p)))

	    (if (null? charts)
		(begin
		  (write-line `(failure after ,(length result) steps))
		  (write-line `(failed at time = ,t seconds))
		  (set! *point-of-failure* p)
		  (error (ode-integrator-error 2))))

	    (let next-chart ((charts charts))
	      (if (null? charts)

		  ;; No more charts: Panic!
		  (begin
		    (write-line `(failure after ,(length result) steps))
		    (write-line `(failed at time = ,t seconds))
		    (set! *point-of-failure* p)
		    (error (ode-integrator-error 1)))

		  ;; Take a step forward in the next chart:
		  (let* ((chart (car charts))
			 (not-in-range
			  (compose not (chart:get-range-test chart)))
			 (v.field (make-local-field chart))
			 (make-field (field-protector v.field)))

		    (let ((new

			   ;; This hack provides an escape mechanism from the
			   ;; local integrator: Check if it tries to access the
			   ;; field at a point outside the current chart.

			   (call-with-current-continuation
			    (lambda (return)
			      (next-step
			       (chart:point->coords p chart)
			       (make-field
				chart
				(list (list not-in-range return))
				'())
			       (lambda () (return #f)))))))

		      (if new
			  (let ((x (integrator:get-new-x new))
				(dt (integrator:get-dt new)))
			    (next-point (chart:coords->point x chart)
					(+ t dt)
					(cons (list t p) result)))
			  (next-chart (cdr charts))))))))
	  result))))
