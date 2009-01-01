;;; This file defines some ways of putting together equations from different
;;; charts.

(declare (usual-integrations))


;;; Just adding linear combinations of old equations to the matrix won't do
;;; anything new.  We also need to eliminate old equations -- How do we do
;;; this?  In the case when two nodes overlap, the choice is simple: Just
;;; eliminate the two original equations, and replace one of the old nodes with
;;; the other.  But what do we do when the node lies *inside* an element?

;;; Let's leave that for future work, and instead implement a method based on
;;; "copying" nodes from other charts to guarantee that nodes overlap exactly
;;; in interactions.  Then this reduces to something like the exact case, and
;;; we know how to combine equations in this case.

;;; By the way, is it important to enforce constraints on boundary nodes?  The
;;; current design makes this impossible to do, but one could probably fix it.

(define (merge-equations domain equations)
  (let ((nodes (manifold:get-nodes domain))
	(count 0)
	(mat #f))

    ;; First, assign IDs to nodes, and create the matrix:

    (write-line '(creating matrix...))

    (let loop ((nodes nodes) (i 0))
      (if (null? nodes)
	  (begin
	    (set! count i)
	    (set! mat (make-sparse-matrix count (+ count 1))))
	  (let ((node (car nodes)))
	    (cond ((node:boundary? node)
		   (node:set-id! node 'boundary-node!)
		   (loop (cdr nodes) i))
		  ((node:get-constraint node)
		   (node:set-id! node 'constrained-node!)
		   (loop (cdr nodes) i))
		  (else
		   (node:set-id! node i)
		   (loop (cdr nodes) (+ i 1)))))))

    ;; Next, start filling in equations while keeping track of constraints:

    (write-line '(copying equations...))

    (let next-eq ((equations equations))
      (if (null? equations)
	  (begin
	    (write-line '(done!))
	    mat)
	  (let* ((eq (car equations))
		 (i (node:get-real-id (equation:get-node eq))))

	    (sparse-matrix-set! mat i count
				(+ (equation:get-constant eq)
				   (sparse-matrix-ref mat i count)))

	    (let next-term ((terms (equation:get-terms eq)))
	      (if (null? terms)
		  (next-eq (cdr equations))
		  (let* ((term (car terms))
			 (j (node:get-real-id (term:get-node term)))
			 (val (term:get-coeff term)))
		    (sparse-matrix-set! mat i j (+ (sparse-matrix-ref mat i j)
						   val))
		    (next-term (cdr terms))))))))))


;;; A slightly different way to merge equations that requires overlaps: Append
;;; linear constraints that force nodal values in overlapping regions to agree
;;; with the interpolated value in the other chart.  Note that if one does this
;;; for *all* nodes, it might "stiffen" the solution over the overlap and force
;;; it to be approximately linear.  Thus, one should try to avoid having too
;;; many constrained nodes, or to somehow reduce the overlap.

(define (append-constraint-equations make-constraints)
  (lambda (domain equations)

    ;; First, set IDs and clear hidden states:

    (write-line '(setting node ids...))

    (let loop ((id 0) (nodes (manifold:get-nodes domain)))
      (if (not (null? nodes))
	  (let ((node (car nodes)))
	    (node:set-constraint! node #f)
	    (if (node:boundary? node)
		(begin
		  (node:set-id! node 'boundary-node!)
		  (loop id (cdr nodes)))
		(begin
		  (node:set-id! node id)
		  (loop (+ id 1) (cdr nodes)))))))

    ;; Next, generate constraints:

    (write-line '(generating constraints...))

    (with-values
	(lambda () (make-constraints domain))

      (lambda (c-count clists)
	(let* ((eq-count (length equations))
	       (m (+ eq-count c-count))
	       (n (+ eq-count 1)))

	  (write-line `(constructing a matrix of dimension (,m ,n)...))

	  (let ((mat (make-sparse-matrix m n)))

	    ;; First, copy the equations:

	    (write-line `(copying ,eq-count equations...))

	    (for-each
	     (lambda (eq)
	       (let ((i (equation:get-id eq)))
		 (sparse-matrix-set!
		  mat i eq-count (equation:get-constant eq))
		 (for-each
		  (lambda (term)
		    (sparse-matrix-set! mat i (term:get-id term)
					(term:get-coeff term)))
		  (equation:get-terms eq))))
	     equations)

	    ;; Next, copy the constraints:

	    (write-line `(copying ,c-count constraints...))

	    (let next-clist ((i eq-count) (clists clists))
	      (if (null? clists)
		  mat
		  (let next-constraint ((clist (car clists)) (i i))
		    (if (null? clist)
			(next-clist i (cdr clists))
			(let ((constraint (car clist)))
			  (sparse-matrix-set!
			   mat i eq-count (equation:get-constant constraint))
			  (for-each
			   (lambda (term)
			     (sparse-matrix-set! mat i (term:get-id term)
						 (term:get-coeff term)))
			   (equation:get-terms constraint))
			  (next-constraint (cdr clist) (+ i 1)))))))))))))


;;; Here's one way to make constraints:

(define (make-ordered-boundary-constraints domain)
  (let* ((charts (manifold:get-finite-atlas domain))
	 (result-1 (charts->constraints charts node:local-boundary?))
	 (result-2 (charts->constraints
		    (reverse charts) node:local-boundary?)))
    (values (+ (car result-1) (car result-2))
	    (append (cadr result-1) (cadr result-2)))))

(define (charts->constraints charts good-node?)
  (let next-chart ((charts charts)
		   (count 0)
		   (clists '()))
    (if (null? charts)

	(list count clists)

	;; Go through each node in the chart and check for constraints:

	(let ((chart (car charts)))
	  (let next-node ((nodes (chart:get-nodes chart))
			  (count count)
			  (clist '()))
	    (if (null? nodes)
		(next-chart (cdr charts) count (cons clist clists))
		(let ((node (car nodes)))
		  (if (and (good-node? node)
			   (not (node:get-constraint node))
			   (not (node:boundary? node)))
		      (let ((eq (make-constraint node (cdr charts))))
			(if eq
			    (next-node (cdr nodes) (+ count 1) (cons eq clist))
			    (next-node (cdr nodes) count clist)))
		      (next-node (cdr nodes) count clist)))))))))

(define (make-constraint node charts)
  (let loop ((charts charts))
    (if (null? charts)
	#f
	(let ((eq (chart:pointwise-constraint node (car charts))))
	  (if eq
	      eq
	      (loop (cdr charts)))))))


;;; A slightly different approach that generates *more* constraints:

(define make-all-ordered-constraints
  (let ((exists? (lambda (node) #t)))
    (lambda (domain)
      (let* ((charts (manifold:get-finite-atlas domain))
	     (result-1 (charts->constraints charts exists?))
	     (result-2 (charts->constraints (reverse charts) exists?)))
	(values (+ (car result-1) (car result-2))
		(append (cadr result-1) (cadr result-2)))))))


;;; Finally, something that generates a lot of constraints:

(define (make-all-constraints domain)
  (let ((constraints
	 (append-map
	  (lambda (pair)
	    (let ((chart-1 (car pair))
		  (chart-2 (cadr pair)))
	      (append (constrain-all-nodes chart-1 chart-2)
		      (constrain-all-nodes chart-2 chart-1))))
	  (pairs (manifold:get-finite-atlas domain)))))
    (values (length constraints) (list constraints))))

(define (constrain-all-nodes chart-1 chart-2)
  (append-map
   (lambda (node)
     (if (node:boundary? node)
	 '()
	 (let ((eq (chart:pointwise-constraint node chart-2)))
	   (if eq
	       (list eq)
	       '()))))
   (chart:get-nodes chart-1)))


;;; Some extensions to charts:

;;; Here's one problem with this approach to constraints, though: Consider the
;;; case when we *do* have a mesh over a simple subset of the plane.  Let's try
;;; to apply this constraint idea to this case: Cut the mesh along some line
;;; formed by the edges, so that we cut the meshed region R into two subregions
;;; R1 and R2.  Let's try to paste R1 and R2 together using constraints.  What
;;; we notice is that when we identify two nodes (by adding a constraint
;;; equation), we are actually requiring that the nodal value satisfies *two*
;;; separate equations, one for R1 and one for R2, instead of satisfying the
;;; *sum* of those two equations.  What this indicates is that this method
;;; actually *requires* overlaps.

;;; Furthermore, note that building constraints by putting in (indeterminate)
;;; Dirichlet boundary conditions along chart boundaries won't work.  Consider
;;; the ideal case, where two charts overlap by exactly their boundary: Since
;;; the Dirichlet problem is well-posed for Laplace's equation, we can put
;;; *any* "boundary value" on this overlap and still get a solution of the
;;; equation over the whole domain.  Clearly, this idea *may* work if one uses
;;; von Neumann conditions rather than Dirichlet conditions, but how to do that
;;; nicely is not clear.  Perhaps using Lagrange multipliers for a constrained
;;; minimization...

(define (chart:pointwise-constraint node chart)

  ;; The coefficients of a linear constraint for some node x should simply be
  ;; the value at p of the basis function centered at x.  This linearity
  ;; depends only on the fact that the solution is approximated by a linear
  ;; combination of basis functions.

  (if (chart:member? (node:get-point node) chart)
      (let* ((x (chart:point->coords (node:get-point node) chart))
	     (element (chart:coords->any-element x chart)))
	(if element
	    (let loop ((nodes (element:get-nodes element))
		       (i 0)
		       (const 0)
		       (terms (list (make-term node -1))))
	      (if (null? nodes)
		  (begin
		    (node:set-constraint! node chart)
		    (make-equation node const terms))
		  (let ((neighbor (car nodes))
			(coeff (evaluate-basis-function
				(element:get-basis-function element i) x)))
		    (if (node:boundary? neighbor)
			(loop (cdr nodes)
			      (+ i 1)
			      (- const (* (node:get-value neighbor) coeff))
			      terms)
			(loop (cdr nodes)
			      (+ i 1)
			      const
			      (cons (make-term neighbor coeff) terms))))))
	    #f))
      #f))
