;;; What do we do when a orbit makes a transition between charts?  This depends
;;; on the representation of points; GJS used an explicit imbedding.  Whatever
;;; the representation, we can use hashing to look up charts given a point (and
;;; this operation should be abstracted anyway).

;;; What about using imbeddings?  What would this do for either ODEs (the paths
;;; should just stay on the manifold) or PDEs?  This is not very general.

;;; Since we can only get global existence under limited circumstances (e.g. a
;;; smooth vector field is complete if the manifold is compact), we might as
;;; well assume that the underlying manifold (or configuration space, since
;;; local existence works for second-order equations as well) is compact, thus
;;; allowing us to assume that the manifold has a finite atlas.  Given that, we
;;; can evolve the integral curve in all charts at the same time, and pick the
;;; "best" solution at each step.  (Or even any solution at all.)

;;; Actually, this idea only requires a locally finite atlas.

(declare (usual-integrations))


;;; Solve ODEs with 4th-order Runge-Kutta:
;;; (Needs automatic step-size control.)


;;; Something to agree on:

(define integrator:package-result cons)
(define integrator:get-new-x cdr)
(define integrator:get-dt car)


;;; Local integrators are assumed to take, in order, the following things:
;;; A STATE vector, a vector field, and FAIL (which is a thunk that does
;;; something in the case of an error).


;;; This procedure turns real vector fields into constructors for local fields;
;;; this exists (mostly) for compatibility issues.

(define (v.field->local-field-maker v.field)
  (lambda (chart)
    (let ((f (chart:get-inverse-map chart)))
      (lambda (x)
	(chart:push-forward (v.field (f x)) chart)))))


;;; RK4 on manifolds.  Still needs QC.

(define (make-rk4-integrator dt)
  (lambda (x v fail)
    (let ((dt/2 (/ dt 2.))
	  (dt/6 (/ dt 6.)))

	    (let* ((F1 (v x))
		   (F2 (v (vector:+ x (vector:* dt/2 F1))))
		   (F3 (v (vector:+ x (vector:* dt/2 F2))))
		   (F4 (v (vector:+ x (vector:* dt F3)))))
	      (integrator:package-result
	       dt
	       (vector:+
		(vector:* dt/6
			  (vector:+ F1
				    (vector:* 2. F2)
				    (vector:* 2. F3)
				    F4))
		x))))))


;;; If ScmUtils is loaded, just use the integrator there (QCRK4 is the default,
;;; I believe):

(define (make-scmutils-integrator dt tol)
  (lambda (x v fail)
    (let* ((w (lambda (x)
		(vector-append (vector 1) (v (vector-tail x 1)))))
	   (result (ode-advancer w (vector-append (vector 0) x) dt tol)))
      (integrator:package-result
       (vector-ref result 0)
       (vector-tail result 1)))))


;;; We now need ODE solvers.  When dealing with ODEs (as opposed to PDEs),
;;; there is one particular problem we need to address: When do we switch
;;; charts, if integrating locally?

;;; Here's an integrator integrates one step at a time.  It understands how and
;;; when to switch between charts.

;;; This procedure uses continuations to handles faults like stepping out of
;;; the chart on a intermediate step, or detecting really bad error conditions.
;;; This simplifies the local integrator.  (Tail recursion is *cool*.)

;;; We should modify this continuation hack to allow QC-RK-4 to punt charts
;;; based on local error analysis.

;;; PRE-CHECKS and POST-CHECKS are predicate-continuation pairs that check for
;;; errors before and after the computation of the vector field.

(define (field-protector v.field)
  (lambda (chart pre-checks post-checks)
    (lambda (x)
      (let loop ((l post-checks))
	(if (null? l)
	    (let ((v (v.field x)))
	      (let loop ((l post-checks))
		(if (null? l)
		    v
		    (let ((pred-cont (car l)))
		      (if ((car pred-cont) v)
			  ((cadr pred-cont) v)
			  (loop (cdr l)))))))
	    (let ((pred-cont (car l)))
	      (if ((car pred-cont) x)
		  ((cadr pred-cont) x)
		  (loop (cdr l)))))))))

;;; As an optimization, it might make sense for vector fields to work directly
;;; with charts.  That is, a vector field is a constructor that, given a chart,
;;; constructs a local vector field on the chart.  This fits in more nicely
;;; with a Lagrangian or Hamiltonian description of mechanics.


;;; There are some problems here:

;;; If we switch charts *before* we step off, then how do we know which chart
;;; to switch to?  The program could sit there switching charts forever unless
;;; we build some memory into this.

;;; On the other hand, if we switch charts *after* we step off, then the
;;; inverse mapping fails.  So what we need is a "compact refinement" of a
;;; covering of charts...

;;; We can always simultaneously evolve the point in several charts.  But is
;;; that too slow?

(define ode-integrator-error
  (let ((errors
	 (vector "*** Warning: Cannot find a good chart.  Using any chart..."
		 "Error: I'm stuck! (Out of charts!) -- V.FIELD->FLOW"
		 "Error: Is this point in the manifold at all?")))
    (lambda (i)
      (vector-ref errors i))))

(define *point-of-failure* #f)

(define (integrator-failure-point)
  *point-of-failure*)


;;; M is the manifold on which we integrate, MAKE-LOCAL-FIELD is a function
;;; that takes a chart and returns a (local) vector field on that chart,
;;; NEXT-STEP is the local integrator, and DISTORTION lets the user rank how
;;; undesirable a particular result is.

(define (v.field->flow M make-local-field next-step distortion)
  (lambda (p t-final . aux)

    ;; Reset the error-reporting mechanism:
    (set! *point-of-failure* #f)

    ;; AUX lets the user specify an initial time (optional).

    (let next-point ((p p)
		     (t (if (not (null? aux)) (car aux) 0.))
		     (result '()))
      (if (<= t t-final)

	  ;; Find all the charts containing this point and try each of them:

	  (let ((charts (manifold:get-local-atlas M p)))

	    ;(write-line '(one more point!))

	    (let next-chart ((charts charts) (min #f) (best #f))
	      (if (null? charts)

		  ;; If we have found a pretty good answer, use it for the next
		  ;; step and save the previous time step.  Otherwise panic and
		  ;; use any chart we can find.  If we can't even find a chart,
		  ;; then the previous step was really bad, too, so the program
		  ;; just dies.

		  (if best

		      (next-point (cadr best)
				  (+ t (car best))
				  (cons (list t p) result))

		      (let ((chart (manifold:find-best-chart M p)))
			(if chart
			    (let* ((v.field (make-local-field chart))
				   (make-field (field-protector v.field))
				   (new (next-step
					 (chart:point->coords p chart)
					 (make-field chart '() '())
					 (lambda () (return #f)))))
			      (newline)
			      (display (ode-integrator-error 0))
			      (newline)
			      (next-point (chart:coords->point
					   (integrator:get-new-x new) chart)
					  (+ t (integrator:get-dt new))
					  (cons (list t p) result)))
			    (begin
			      (write-line
			       `(failure after ,(length result) steps))
			      (write-line `(failed at time = ,t seconds))
			      (set! *point-of-failure* p)
			      (error (ode-integrator-error 1))))))

		  ;; Take a step forward in the next chart:

		  (let* ((chart (car charts))
			 (not-in-range
			  (compose not (chart:get-range-test chart)))
			 (v.field (make-local-field chart))
			 (make-field (field-protector v.field)))

		    ;(write-line '(one more chart!))

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

		      ;; If nothing went wrong, check if the current chart does
		      ;; better, and keep the result if it does.  Otherwise
		      ;; keep the old results.

		      (if new
			  (let ((x (integrator:get-new-x new)))
			    (if (chart:in-range? x chart)
				(let* ((q (chart:coords->point x chart))
				       (dt (integrator:get-dt new))
				       (e (distortion
					   chart
					   (make-tangent chart q
							 (v.field x)))))
				  ;(write-line `(t = ,t e = ,e))
				  (if min
				      (if (< e min)
					  (next-chart (cdr charts)
						      e
						      (list dt q))
					  (next-chart (cdr charts) min best))
				      (next-chart (cdr charts) e (list dt q))))
				(next-chart (cdr charts) min best)))
			  (next-chart (cdr charts) min best)))))))
	  result))))


;;; Local integrators, on the other hand, are more general.  But they have the
;;; following problems:
;;;
;;; 1. They are less efficient.
;;;
;;; 2. When to switch charts?
;;;
;;; 3. Distortion of the vector field may produce locally bad solutions.

;;; Some of these problems (namely 2) can be solved by requiring that manifolds
;;; have charts (U, V, f) where the closure of U and V are both compact, and f
;;; extends to a diffeomorphism between those compact sets.  This lets us
;;; switch charts without having a measure of how badly the vector field's
;;; doing.  (Namely, switch when we step out of a chart!)

;;; Using the imbedding to integrate: It's probably faster, in many cases more
;;; intuitive, and avoids the issue of switching between charts.  The problems
;;; are:
;;;
;;; 1. It's not as useful in abstract manifolds.  We don't know if all
;;;    manifolds important for applications will be represented by imbeddings).
;;;
;;; 2. The trajectories may not stay on the manifold.

;;; Some of these problems can be resolved by having a uniform way to attach
;;; special structures to charts.  It's ugly, but it'll be needed for solving
;;; PDEs.
