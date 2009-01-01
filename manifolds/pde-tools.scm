;;; This file provides the interface between the manifold PDE code and the FEM
;;; toolkit.  These should be actual examples of constructors for PDE charts
;;; that use the FEM stuff we already have.  This part can be as specific as is
;;; necessary.

;;; If we didn't care about good triangulations, then one might think that this
;;; recursive algorithm works for convex sets: Sort the nodes along one axis,
;;; and form the n-simplex with the highest n nodes.  Take the vertices not in
;;; contact with the rest of the solid out, and apply recursively.  This is an
;;; O(n log n) algorithm.  Well, it actually doesn't seem to work, not without
;;; some tweaking.  In any case, there is something better:

;;; See http://www.geom.umn.edu/software/qhull/ for general triangulation
;;; algorithms; it's probably better than the Guibas/Stolfi algorithm we're
;;; using.  (These algorithms usually work by requiring convexity, which we
;;; will also do.)

(declare (usual-integrations))


;;; Local tesselations depend on the fact that domains of charts are often
;;; simple shapes (like squares and discs).  Focus on 2-D case for now.

(define make-rectangular-vertices
  (let ((border-frac 1e-3))
    (lambda (chart x-low x-high y-low y-high vcount hcount)

      ;; Stay away from the borders so that all the nodes lie in the chart.

      (let ((x-border (* (- x-high x-low) border-frac))
	    (y-border (* (- y-high y-low) border-frac)))
	(set! x-high (- x-high x-border))
	(set! x-low (+ x-low x-border))
	(set! y-high (- y-high y-border))
	(set! y-low (+ y-low y-border))

	;; Is this a boundary chart?  If so, fix it so that boundary nodes are
	;; really on the boundary:

	(if (boundary-chart? chart)
	    (let ((level (chart:get-boundary-index chart)))
	      (if (= (chart:get-boundary-index chart) 0)
		  (if (almost-equal? level (- x-low x-border))
		      (set! x-low (- x-low x-border))
		      (set! x-high (+ x-high x-border)))
		  (if (almost-equal? level (- y-low y-border))
		      (set! y-low (- y-low y-border))
		      (set! y-high (+ y-high y-border)))))))

      (let* ((hcount-1 (- hcount 1))
	     (vcount-1 (- vcount 1))
	     (dx (/ (- x-high x-low) hcount-1))
	     (dy (/ (- y-high y-low) vcount-1)))

	(let next-row ((i 0) (nodes '()))
	  (if (< i hcount)
	      (let ((x (+ x-low (* i dx))))
		(let next-col ((j 0) (nodes nodes))
		  (if (< j vcount)
		      (let* ((y (+ y-low (* j dy)))
			     (new-node (make-node (vector x y) chart)))
			(if (or (= i 0) (= i hcount-1) (= j 0) (= j vcount-1))
			    (node:set-local-boundary! new-node #t))
			(next-col (+ j 1) (cons new-node nodes)))
		      (next-row (+ i 1) nodes))))
	      nodes))))))

(define make-circular-vertices
  (let ((border-frac 1e-3))
    (lambda (chart x y radius radial-count angular-count)

      ;; Stay away from the borders so that all the nodes lie in the chart.  By
      ;; our convention, this can't be a boundary chart, so no boundary nodes.

      (set! radius (* (- 1 border-frac) radius))

      (let* ((radial-count-1 (- radial-count 1))
	     (dr (/ radius radial-count-1))
	     (dt (/ (* 2 pi) angular-count)))

	(let next-ray ((i 1) (nodes '()))
	  (if (< i radial-count)
	      (let ((r (* i dr)))
		(let next-point ((j 0) (nodes nodes))
		  (if (< j angular-count)
		      (let* ((t (* j dt))
			     (new-node (make-node (vector (+ (* r (cos t)) x)
							  (+ (* r (sin t)) y))
						  chart)))
			(if (= i radial-count-1)
			    (node:set-local-boundary! new-node #t))
			(next-point (+ j 1) (cons new-node nodes)))
		      (next-ray (+ i 1) nodes))))
	      (cons (make-node (vector x y) chart) nodes)))))))


;;; Let's use the range structures for meshing.  Better back off from the edge
;;; a bit to a compact subset, though, to avoid trouble with charts that cover
;;; almost all of the manifold (such as spherical coordinates on the sphere).

(define (make-vertices chart . args)
  (let ((dim (chart:dimension chart)))
    (if (= dim 2)
	(let ((nodes
	       (cond ((chart:cell-range? chart)
		      (apply make-rectangular-vertices
			     (append
			      (cons chart
				    (append-map
				     (lambda (interval)
				       (list (interval:inf interval)
					     (interval:sup interval)))
				     (cell-range:get-interval-list chart)))
			      (cdr (assq 'rectangular args)))))

		     ((chart:spherical-range? chart)
		      (let ((center (spherical-range:get-center chart)))
			(apply make-circular-vertices
			       (append
				(list chart
				      (vector-ref center 0)
				      (vector-ref center 1)
				      (spherical-range:get-radius chart))
				(cdr (assq 'spherical args))))))
		     (else (error "Don't know how to mesh this chart!")))))
	  nodes)
	(write-line '(can only handle planar regions for now!)))))


;;; A discretization routine.  This is our interface to local FEM.

(define (fem-discretize chart source)
  (let* ((nodes (list->vector (chart:get-nodes chart)))
	 (sparse (assemble-equations source nodes))
	 (n (sparse-matrix-row-count sparse))
	 (m (vector-length nodes))
	 (index-map (make-vector n #f)))

    (write-line `(,n equations generated for ,m nodes.))

    ;; Need to be able to translate row indices into node IDs.

    (let loop ((i 0) (j 0))
      (if (< i m)
	  (let ((node (vector-ref nodes i)))
	    (if (node:boundary? node)
		(loop (+ i 1) j)
		(begin
		  (vector-set! index-map j i)
		  (loop (+ i 1) (+ j 1)))))))

    (let loop ((i 0) (result '()))
      (if (< i n)
	  (let ((node (vector-ref nodes (vector-ref index-map i))))
	    (let next-term ((row (sparse-matrix-get-row sparse i))
			    (terms '())
			    (const 0))
	      (if (null? row)
		  (loop (+ i 1) (cons (make-equation node const terms) result))
		  (let* ((pair (car row))
			 (index (car pair))
			 (coeff (cadr pair)))
		    (if (= index n)
			(next-term (cdr row) terms coeff)
			(let ((node (vector-ref
				     nodes (vector-ref index-map index))))
			  (next-term (cdr row)
				     (cons (make-term node coeff) terms)
				     const)))))))
	  result))))


;;; Usually, we won't need more than just the vertices:

(define (make-no-extra-nodes complex)
  (map (lambda (face) '()) (complex->faces complex)))


;;; However, it's always nice to have mid-edge nodes:

(define make-mid-edge-nodes
  (let ((edge-index
	 (lambda (edge)
	   (apply symmetric->vector-index (map node:get-id edge)))))
    (lambda (complex)
      (let* ((nodes (complex->vertices complex))
	     (edges (complex->edges complex))
	     (n (choose (+ (length nodes) 2) 2))
	     (big-v (make-vector n #f)))

	;; First, assign each node a unique ID:

	(let loop ((i 0) (nodes nodes))
	  (if (not (null? nodes))
	      (begin
		(node:set-id! (car nodes) i)
		(loop (+ i 1) (cdr nodes)))))

	;; Next, construct mid-edge nodes for each edge and save them:

	(for-each
	 (lambda (edge)
	   (let* ((org (car edge))
		  (dest (cadr edge))
		  (node (make-node
			 (vector:* 1/2
				   (apply vector:+ (map node:get-coords edge)))
			 (node:get-chart org))))

	     (if (and (node:local-boundary? org)
		      (node:local-boundary? dest))
		 (node:set-local-boundary! node #t))

	     (vector-set! big-v (edge-index edge) node)))

	 edges)

	;; Finally, assign edge nodes to each face of the complex:

	(map
	 (lambda (face)
	   (append-map
	    (lambda (pair)
	      (let ((node (vector-ref big-v (edge-index pair))))
		(if node
		    (list node)
		    '())))
	    (pairs face)))
	 (complex->faces complex))))))

(define (complex->edges complex)
  (cadr (reverse complex)))

(define (complex->vertices complex)
  (car (reverse complex)))

(define (complex->faces complex)
  (car complex))


;;; Mesh-generation stuff: Here's an interface to Delaunay triangulation
;;; routines in 2-D.

(define (planar-triangulate nodes)
  (reverse (cons nodes (delaunay-triangulation (list->vector nodes)))))


;;; Count the number of boundary nodes in a VECTOR of nodes:

(define (number-of-boundary-nodes nodes)
  (let ((n (vector-length nodes)))
    (let loop ((i 0) (count 0))
      (if (< i n)
	  (if (node:boundary? (vector-ref nodes i))
	      (loop (+ i 1) (+ count 1))
	      (loop (+ i 1) count))
	  count))))
