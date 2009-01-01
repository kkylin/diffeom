;;; This program relies on the existence of local triangulation algorithms for
;;; arbitrary dimensions.  That problem does appear to be solved: See Barber,
;;; Dobkin, Huhdanpaa, "The Quickhull Algorithm for Convex Hulls."  The paper,
;;; as well as the software itself, are available at the University of
;;; Minnesota's Geometry Center, "http://www.geom.umn.edu/software/qhull/".

;;; See pde-main.scm.old for lots and lots of comments and design notes.  (It's
;;; a bit incoherent, so I've deleted them from this version.)

;;; NOTE:

;;; Yet another possible approach is to look through the differential topology
;;; literature and see if there exists a (constructive) proof that manifolds
;;; are triangulable or CW complexes.  One can probably show that manifolds are
;;; CW complexes by Morse theory; is this computable?  Does "computational
;;; Morse theory" exist?

;;; If we require the manifolds to be CW complexes or simplicial compexes to
;;; begin with, all these problems would be solved.

(declare (usual-integrations))
(load "pde-mergers")

;;; Given a domain with constructed elements, a source function, and a boundary
;;; value function, produce the appropriate discretized equation.  The nodes
;;; are left with indices that specify their corresponding row in the matrix.

(define (pde:equation-maker merge-equations)
  (lambda (domain source boundary-value . extra-args)

    ;; EXTRA-ARGS gives us finer control over the discretization.

    ;; DOMAIN should be a manifold that already has PDE structures constructed.
    ;; Hence, it contains information about the operator (through the elements
    ;; in its discretized charts).

    ;; BOUNDARY-VALUE is irrelevant for domains without boundary.  Just specify
    ;; anything (but do put in something).

    (let* ((M domain)
	   (charts (manifold:get-finite-atlas M))
	   (nodes (list->vector (append-map chart:get-nodes charts)))
	   (ncount (vector-length nodes)))

      ;; CHART:DISCRETIZE-PDE should return a list of linear equations.  First,
      ;; set the boundary values:

      (write-line `(,ncount nodes generated...))
      (write-line '(setting boundary values...))

      (do ((i 0 (+ i 1)))
	  ((>= i ncount))
	(let ((node (vector-ref nodes i)))
	  (if (node:boundary? node)
	      (node:set-value! node (boundary-value node)))))

      ;; Next, compute the local equation systems:

      (write-line `(computing ,(length charts) local systems of equations...))

      (let ((equations (append-map
			(lambda (chart)
			  (chart:discretize-pde chart source extra-args))
			charts)))

	;; Compute constraints:

	(write-line '(merging local equations...))
	(show-time
	 (lambda ()
	   (merge-equations domain equations)))))))


;;; This procedure creates a constructor that, given a manifold M, creates the
;;; data structures necessary for solving PDEs.  The arguments must agree with
;;; the manifold on a contract that lets procedures obtain chart information.

(define (pde:domain-maker generate-node-lists process-complex)
  (lambda (M
	   make-vertices
	   make-extra-nodes
	   tesselate
	   . argl)

    ;; First, make the bounding nodes of the convex domain, and then
    ;; triangulate and make the extra nodes:

    (let ((atlas (manifold:get-finite-atlas M)))

      (if (not atlas)
	  (error "Error: Can only do FEM with finite atlases."))

      (write-line '(tesselating domain...))

      ;; Do something more complicated here to reduce the overlap:

      (let loop ((charts atlas)
		 (node-lists (generate-node-lists make-vertices atlas argl)))
	(if (not (null? charts))

	    ;; TESSELATE should return a list of lists, where each list
	    ;; contains the elemental faces of a given dimension (in some given
	    ;; polytope).  In the planar case, this reverses the convention in
	    ;; fem.scm: The list should be sorted by dimension in *descending*
	    ;; order.

	    (let* ((chart (car charts))
		   (nodes (car node-lists))
		   (complex (process-complex (tesselate nodes) (cdr charts)))
		   (extra-nodes (make-extra-nodes complex)))

	      ;; By default, use FEM-DISCRETIZE.  Can replace with others.

	      (make-pde-chart chart extra-nodes fem-discretize complex)
	      (loop (cdr charts) (cdr node-lists)))))

      ;; Construct elements.  We don't need to explicitly mark boundaries
      ;; because manifolds should already have such structures defined.

      (lambda (operator make-integrator make-basis-function)
	(let ((element-maker (pde:element-maker operator
						make-integrator
						make-basis-function)))

	  (write-line '(constructing elements...))

	  (for-each

	   (lambda (chart)

	     ;; Construct the elements:

	     (write-line
	      `(making ,(length (complex->faces (chart:get-complex chart)))
		       elements...))

	     (let* ((make-element (element-maker chart))
		    (new-elements (show-time
				   (lambda ()
				     (map make-element
					  (complex->faces
					   (chart:get-complex chart))
					  (chart:get-extra-nodes chart))))))
	       (chart:set-elements! chart new-elements)))

	   atlas))))))


;;; This procedure uses the implicit ordering of the charts to remove extra
;;; nodes in overlaps and to duplicate enough nodes so that the meshes can be
;;; "glued" together.

(define (make-nodes-for-each-chart make-nodes charts extra-args)
  (map (lambda (chart) (apply make-nodes (cons chart extra-args))) charts))

(define (generate-node-lists make-nodes charts argl)

  ;; Generate a list of nodes for each chart, then loop over the charts.  Note
  ;; that the earlier a chart is in the list, the less likely its nodes are to
  ;; survive.

  (let next-chart ((charts charts)
		   (lists (make-nodes-for-each-chart make-nodes charts argl))
		   (result '())
		   (reversed '())
		   (count 0))
    (if (null? charts)
	(copy-overlap-nodes count result reversed)
	(next-chart (cdr charts)
		    (cdr lists)
		    (cons (remove-overlap-nodes (car lists) (cdr charts))
			  result)
		    (cons (car charts) reversed)
		    (+ count 1)))))

(define (copy-between-node-lists make-nodes charts argl)

  ;; Same as GENERATE-NODE-LISTS, but doesn't call REMOVE-OVERLAP-NODES.

  (let ((node-lists (make-nodes-for-each-chart make-nodes charts argl)))
    (copy-overlap-nodes (length charts)
			(reverse node-lists)
			(reverse charts))))

;;; Take out all nodes in NODES that belong to any of the charts in CHARTS.

(define (remove-overlap-nodes nodes charts)
  (let next-node ((nodes nodes) (result '()))
    (if (null? nodes)
	result
	(let* ((node (car nodes))
	       (p (node:get-point node)))
	  (let next-chart ((charts charts))
	    (if (null? charts)
		(next-node (cdr nodes) (cons node result))
		(if (chart:member? p (car charts))
		    (next-node (cdr nodes) result)
		    (next-chart (cdr charts)))))))))

;;; For each node list in LISTS, take each node and see if it's in one of the
;;; charts that come after the node's own chart in list-order.  If so, make a
;;; copy of that node and put it in the corresponding chart.  Note that the
;;; order of node lists is reversed.

(define (copy-overlap-nodes count lists charts)
  (let ((v (make-vector count '())))
    (let next-list ((lists lists) (charts charts) (i 0) (result '()))
      (if (null? lists)
	  result
	  (let next-node ((nodes (car lists)))
	    (if (null? nodes)
		(next-list (cdr lists) (cdr charts) (+ i 1)
			   (cons (append (vector-ref v i) (car lists)) result))
		(let ((node (car nodes)))
		  (if (or (node:local-boundary? node)
			  (node:boundary? node))
		      (let ((p (node:get-point node)))
			(let next-chart ((charts (cdr charts))
					 (j (+ i 1))
					 (l (cdr lists)))
			  (if (null? charts)
			      (next-node (cdr nodes))
			      (let ((chart (car charts)))
				(if (chart:member? p chart)
				    (let ((other (close-node p (car l))))
				      (if other
					  (node:set-constraint! other node)
					  (vector-set! v j
						       (cons
							(node:copy node chart)
							(vector-ref v j))))))
				(next-chart (cdr charts) (+ j 1) (cdr l))))))
		      (next-node (cdr nodes))))))))))

(define close-node
  (let* ((close-enuf? (make-comparator .01))
	 (too-close? (lambda (p q)
		       (close-enuf? (vector:distance p q) 0))))
    (lambda (p l)
      (let loop ((l l))
	(if (null? l)
	    #f
	    (if (too-close? p (node:get-point (car l)))
		(car l)
		(loop (cdr l))))))))

;;; After filtering out nodes, local boundary information becomes useless...

(define (exact-overlap complex charts)
  (kill-extra-nodes complex charts)
  (resurrect-only-connected-nodes complex charts)
  (keep-only-live-nodes complex charts))

(define (remove-overlap complex charts)
  (kill-extra-nodes complex charts)
  (resurrect-some-nodes complex charts)
  (keep-only-live-nodes complex charts))

(define (reduce-overlap complex charts)
  (kill-extra-nodes complex charts)
  (resurrect-some-nodes complex charts)
  (resurrect-some-nodes complex charts)
  (keep-only-live-nodes complex charts))

(define (absolutely-no-overlap complex charts)
  (kill-extra-nodes complex charts)
  (keep-only-live-nodes complex charts))

(define (extended-overlap complex charts)
  (extend-local-boundary complex charts)
  complex)

(define (extend-local-boundary complex charts)

  (write-line '(extending local boundary...))

  (let loop ((edges (complex->edges complex)) (keep '()))
    (if (null? edges)
	(for-each
	 (lambda (node)
	   (node:set-local-boundary! node #t))
	 keep)
	(let ((n1 (caar edges))
	      (n2 (cadar edges)))
	  (if (node:local-boundary? n1)
	      (if (node:local-boundary? n2)
		  (loop (cdr edges) keep)
		  (loop (cdr edges) (cons n2 keep)))
	      (if (node:local-boundary? n2)
		  (loop (cdr edges) (cons n1 keep))
		  (loop (cdr edges) keep)))))))

(define (do-nothing-to-complex complex charts)
  complex)

(define (kill-extra-nodes complex charts)

  ;; Figure out which nodes to keep by looking at the overlaps:

  (write-line `(processing ,(length (complex->vertices complex)) nodes...))

  (let next-node ((nodes (complex->vertices complex)))
    (if (not (null? nodes))
	(let ((node (car nodes)))
	  (let ((p (node:get-point node)))
	    (let next-chart ((charts charts))
	      (if (null? charts)
		  (next-node (cdr nodes))
		  (let ((chart (car charts)))
		    (if (chart:member? p chart)
			(let ((node (car nodes)))
			  (node:kill! node)
			  (node:set-local-boundary! node #f)
			  (next-node (cdr nodes)))
			(next-chart (cdr charts)))))))))))

(define (resurrect-some-nodes complex charts)

  ;; We actually need to keep more nodes than this, because we need *some*
  ;; overlap, though not too much.  We also need more sophisticated ways of
  ;; checking whether a node should be kept.

  (write-line '(figuring out overlaps...))

  (let loop ((faces (complex->faces complex)) (keep '()))
    (if (null? faces)
	(for-each
	 (lambda (face)
	   (for-each
	    (lambda (node)
	      (if (not (node:active? node))
		  (begin
		    (node:set-local-boundary! node #t)
		    (node:resurrect! node))))
	    face))
	 keep)
	(if (save-face? (car faces) charts)
	    (loop (cdr faces) (cons (car faces) keep))
	    (loop (cdr faces) keep)))))

(define (resurrect-only-connected-nodes complex charts)

  ;; Only keep nodes that are connected to live ones:

  (write-line '(figuring out overlaps...))

  (let loop ((faces (complex->faces complex)) (keep '()))
    (if (null? faces)
	(for-each
	 (lambda (face)
	   (for-each
	    (lambda (node)
	      (if (not (node:active? node))
		  (begin
		    (node:set-local-boundary! node #t)
		    (node:resurrect! node))))
	    face))
	 keep)
	(if (at-least-one-live-node? (car faces) charts)
	    (loop (cdr faces) (cons (car faces) keep))
	    (loop (cdr faces) keep)))))

(define (keep-only-live-nodes complex charts)

  ;; Figure out which faces/edges/etc. to keep:

  (write-line '(processing complex...))

  (let loop ((complex complex) (result '()))
    (if (null? complex)
	(reverse result)
	(let inner-loop ((faces (car complex)) (okay-faces '()))
	  (if (null? faces)
	      (loop (cdr complex) (cons okay-faces result))
	      (let* ((face (car faces))
		     (list? (list? face)))
		(if (or (and list? (not (memq #f (map node:active? face))))
			(and (not list?) (node:active? face)))
		    (inner-loop (cdr faces) (cons face okay-faces))
		    (inner-loop (cdr faces) okay-faces))))))))


;;; We need to check if the particular element covers a region of the manifold
;;; that isn't covered by another chart.  It is hard to do this in general, but
;;; we can use a probabilistic algorithm and take advantage of the fact that
;;; elements are convex:

(define save-face?

  ;; KLUGE-FACTOR defines how many random points to try.  It should scale up
  ;; with the dimension/size of the element, but for now it's constant.

  (let ((kluge-factor 40))
    (lambda (face charts)

      (or (at-least-one-live-node? face charts)

	  ;; Go through a complicated (and probabilistic) test:

	  (let ((vertices (map node:get-coords face))
		(chart (node:get-chart (car face))))
	    (let ((m (length face))
		  (n (chart:dimension chart)))
	      (let loop ((k 0))
		(if (< k kluge-factor)

		    ;; Contruct a random point in the element:

		    (let ((v (make-random-probability-vector m))
			  (w (make-vector n 0)))

		      (do ((i 0 (+ i 1))
			   (vertices vertices (cdr vertices)))
			  ((>= i m))
			(let ((coeff (vector-ref v i))
			      (x (car vertices)))
			  (do ((j 0 (+ j 1)))
			      ((>= j n))
			    (vector-set! w j (+ (vector-ref w j)
						(* coeff (vector-ref x j)))))))

		      ;; Now test it:

		      (let ((p (chart:coords->point w chart)))
			(let next-chart ((charts charts))
			  (if (null? charts)
			      #t
			      (if (chart:member? p (car charts))
				  (loop (+ k 1))
				  (next-chart (cdr charts)))))))
		    #f))))))))

(define (make-random-probability-vector n)
  (let ((v (make-vector n)))
    (let loop ((i 0) (sum 0))
      (if (< i n)
	  (let ((val (random 1.)))
	    (vector-set! v i val)
	    (loop (+ i 1) (+ val sum)))
	  (do ((j 0 (+ j 1)))
	      ((>= j n) v)
	    (vector-set! v j (/ (vector-ref v j) sum)))))))


;;; Here is a simpler variant:

(define (at-least-one-live-node? face charts)
  (memq #t (map node:active? face)))


;;; A useful procedure that gets all the nodes out of the domain:

(define (manifold:get-nodes domain)
  (append-map chart:get-nodes (manifold:get-finite-atlas domain)))
