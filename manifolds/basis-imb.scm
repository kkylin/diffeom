;;; The basis functions defined here are much like polynomial basis functions,
;;; only they exist directly on the imbedding representation of a manifold,
;;; instead of on the chart.  Many of the procedures in 2d-poly-basis.scm are
;;; called here.

(declare (usual-integrations))


;;; Interface to the manifold code:

(define (pde:make-imbedded-poly-basis-function nodes i)
  (let ((basis (make-imbedded-basis-function
		nodes i (node:get-chart (car nodes)))))
    (node:add-basis-function (list-ref nodes i) basis)
    basis))

(define (operator:imbedded-poly-op left-op right-op combine)
  (lambda (chart nodes)
    (lambda (f g)
      (combine (left-op f) (right-op g)))))


;;; Basic constructor:

(define (vector->imbedded v chart)
  (package-basis-function-methods
   '2d-imbeded-basis-function
   v
   (imbedded->function v chart)
   (make-2d-poly-adder v)
   (make-2d-poly-subtractor v)
   (make-2d-poly-multiplier v)
   (make-2d-poly-scalar-multiplier v)))

(define (make-imbedded-basis-function nodes center chart)
  (let* ((n (length nodes))
	 (vals (make-vector n))
	 (points (make-vector n)))

    (let loop ((nodes nodes) (i 0))
      (if (null? nodes)
	  (vector->imbedded (poly:point-value->coeff vals points) chart)
	  (let ((node (car nodes)))
	    (if (= i center)
		(vector-set! vals i 1)
		(vector-set! vals i 0))
	    (vector-set! points i (node:get-point node))
	    (loop (cdr nodes) (+ i 1)))))))


;;; A slightly different kind of constructor:

(define (function->imbedded f nodes)
  (let* ((n (length nodes))
	 (vals (make-vector n))
	 (points (make-vector n)))

    (let loop ((i 0) (nodes nodes))
      (if (null? nodes)
	  (vector->poly (poly:point-value->coeff vals points))
	  (let ((node (car nodes)))
	    (vector-set! points i (node:get-point node))
	    (vector-set! vals i (f node))
	    (loop (+ i 1) (cdr nodes)))))))


;;; And its inverse:

(define (imbedded->function f chart)
  (lambda (x)
    (imbedded:evaluate f x chart)))

(define (imbedded:evaluate f x chart)
  (vector-first (poly:coeff->point-value
		 f (vector (chart:coords->point x chart)))))


;;; The truly messy stuff: Integrals!  This needs to run a lot faster.  What
;;; about doing away with the coordinate transformations?

(define (make-triangular-imbedded-integrator vertex-nodes)

  ;; We assume that there are three vertex nodes, and that the triangle they
  ;; form is the boundary of the element:

  (if (not (= (length vertex-nodes) 3))
      (error (string-append "Error: Elements must have three vertex nodes."
			    " -- MAKE-TRIANGULAR-IMBEDDED-INTEGRATOR")))

  (let ((p1 (car vertex-nodes))
	(p2 (cadr vertex-nodes))
	(p3 (caddr vertex-nodes)))

    ;; Find the absolute value of the Jacobian of the affine transformation
    ;; mapping the reference triangle {(0,0),(1,0),(0,1)} to this triangle.

    (let* ((A (list->matrix
	       2 2
	       (list
		(- (node:get-real-x p2) (node:get-real-x p1))
		(- (node:get-real-x p3) (node:get-real-x p1))
		(- (node:get-real-y p2) (node:get-real-y p1))
		(- (node:get-real-y p3) (node:get-real-y p1)))))
	   (b (node:get-point p1))
	   (jacobian (abs (det A))))

      (define (integrate f . rest)
	(let* ((f (apply basis:* (cons f rest)))
	       (degree (poly:degree f))
	       (reference (poly:make-sample-points degree))
	       (n (choose (+ degree 2) 2))
	       (real (make-vector n)))

	  (do ((i 0 (+ i 1)))
	      ((>= i n))
	    (vector-set! real i
			 (apply-affine-transformation
			  A b (vector-ref reference i))))

	  (* jacobian
	     (inner-product
	      (poly:point-value->coeff
	       (poly:coeff->point-value (basis:get-rep f) real) reference)
	      (make-reference-integrals degree)))))

      integrate)))
