;;; This file uses the Delaunay triangulation code to construct two-dimensional
;;; domains of solution.  It defines constructors for different types of
;;; domains of solution, and also defines ways of setting different kinds of
;;; boundary conditions.

(declare (usual-integrations))


;;; A generic template for making domains.  Note that this is limited to planar
;;; regions.  For higher dimensions, we need to handle higher-dimensional
;;; simplices, not just faces and edges...

(define (domain-maker make-vertex-nodes
		      make-edge-nodes
		      make-interior-nodes
		      tesselate
		      make-element
		      make-boundary)

  ;; MAKE-VERTEX-NODES is given the arguments to MAKE-DOMAIN, and should return
  ;; a vector of (vertex) nodes.

  ;; MAKE-EDGE-NODES should return a list of edge nodes on the edge defined by
  ;; the two given nodes.

  ;; MAKE-INTERIOR-NODES should make interior nodes in the element defined by a
  ;; list of given nodes.

  ;; TESSELATE takes a vector of nodes and returns a pair containing a list of
  ;; edges and a list of faces, in that order.  (Edges are defined by pairs of
  ;; nodes, while faces are defined by triplets of nodes.  Just like abstract
  ;; simplicial complexes...)

  ;; MAKE-BOUNDARY takes a vector of nodes and sets the appropriate ones to be
  ;; boundary nodes.

  ;; See fem.scm for the definition of MAKE-ELEMENT.

  (define (make-domain . args)

    ;; Make vertex nodes:

    (write-line '(making nodes...))

    (let* ((vertex-nodes (apply make-vertex-nodes args))
	   (vertex-count (vector-length vertex-nodes)))

      ;; Tesselate the vertex nodes:

      (write-line '(tesselating nodes...))

      (let* ((complex (tesselate vertex-nodes))
	     (edges (car complex))
	     (faces (cadr complex)))

	;; Record some debugging information:

	(set! *debugging-info* edges)

	;; Create edge nodes:

	(write-line '(creating edge nodes...))

	(let ((edge-nodes (make-vector (choose (+ vertex-count 1) 2) '()))
	      (edge-index symmetric->vector-index))

	  (do ((i 0 (+ i 1)))
	      ((>= i vertex-count))
	    (node:set-id! (vector-ref vertex-nodes i) i))

	  (for-each

	   (lambda (pair)
	     (vector-set! edge-nodes
			  (apply edge-index (map node:get-id pair))
			  (apply make-edge-nodes pair)))

	   edges)

	  ;; Create interior nodes:

	  (write-line '(creating interior nodes and making elements...))

	  (let loop ((faces faces) (interior-nodes '()))
	    (if (not (null? faces))
		(let* ((face (car faces))
		       (elist (append-map
			       (lambda (pair)
				 (vector-ref edge-nodes
					     (apply edge-index pair)))
			       (pairs (map node:get-id face))))
		       (ilist (apply make-interior-nodes (append face elist))))

		  (make-element face (append elist ilist))
		  (loop (cdr faces) (append ilist interior-nodes)))

		(begin

		  ;; We now need to combine the three lists of nodes into one
		  ;; big list, and to figure out the boundary:

		  (write-line '(cleaning up...))

		  ;; First, count the number of edge nodes and create a vector
		  ;; to store the node:

		  (let ((edge-count (vector-length edge-nodes)))
		    (let loop ((count 0) (i 0))
		      (if (< i edge-count)
			  (loop (+ count (length (vector-ref edge-nodes i)))
				(+ i 1))
			  (let* ((icount (length interior-nodes))
				 (nodes (make-vector
					 (+ vertex-count count icount))))

			    ;; Report data:

			    (write-line `(,count edge nodes))
			    (write-line `(,icount interior nodes))
			    (write-line `(,vertex-count vertices))

			    ;; Copy the vertex nodes:

			    (write-line '(copying vertices...))

			    (do ((i 0 (+ i 1)))
				((>= i vertex-count))
			      (vector-set! nodes i
					   (vector-ref vertex-nodes i)))

			    ;; Copy the edge nodes:

			    (write-line '(copying edge nodes...))

			    (let loop ((i 0) (j vertex-count))
			      (if (< i edge-count)
				  (loop (+ i 1)
					(let loop ((l (vector-ref edge-nodes
								  i))
						   (j j))
					  (if (null? l)
					      j
					      (begin
						(vector-set! nodes j (car l))
						(loop (cdr l) (+ j 1))))))))

			    ;; Copy the interior nodes:

			    (write-line '(copying interior nodes...))

			    (let loop ((i (+ vertex-count count))
				       (l interior-nodes))
			      (if (null? l)
				  (begin

				    ;; Make boundary nodes:

				    (write-line '(setting boundary...))
				    (make-boundary nodes)

				    ;; Sort and return nodes:

				    (write-line '(sorting nodes...))
				    (sort! nodes lexicographic<))

				  (begin
				    (vector-set! nodes i (car l))
				    (loop (+ i 1) (cdr l))))))))))))))))

  make-domain)


;;; Some useful routines:

(define (make-no-edge-nodes v0 v1) '())

(define (make-no-interior-nodes . args) '())

(define (predicate->make-boundary boundary?)
  (lambda (nodes)
    (let ((n (vector-length nodes)))
      (do ((i 0 (+ i 1)))
	  ((>= i n))
	(let ((node (vector-ref nodes i)))
	  (node:set-boundary! node (boundary? node)))))))

(define (make-midpoint-node v0 v1)
  (list (make-node (/ (+ (node:get-x v0) (node:get-x v1)) 2.)
		   (/ (+ (node:get-y v0) (node:get-y v1)) 2.))))

(define (do-nothing-to-nodes nodes)
  'done)


;;; Various procedures to help build domains:

;;; Circular domains:

(define (make-circular-domain-vertices angular-count radial-count)
  (let* ((rcount+1 (+ radial-count 1))
	 (count (+ (* rcount+1 angular-count) 1))
	 (nodes (make-vector count))
	 (dt (/ (* 2 3.141592653589793) angular-count))
	 (dr (/ .5 rcount+1)))

    (vector-set! nodes 0 (make-node .5 .5))

    (let next-ray ((i 0) (count 1))
      (if (< i angular-count)
	  (let ((t (* i dt)))
	    (let next-node ((j 1) (count count))
	      (if (<= j rcount+1)
		  (let ((r (* j dr)))
		    (vector-set! nodes count
				 (make-node (+ .5 (* r (cos t)))
					    (+ .5 (* r (sin t)))))
		    (next-node (+ j 1) (+ count 1)))
		  (next-ray (+ i 1) count))))
	  nodes))))

(define (circular-boundary? node)
  (let ((x (node:get-x node))
	(y (node:get-y node)))
    (almost-zero? (- .5 (sqrt (+ (square (- x .5)) (square (- y .5))))))))


;;; Square domains:

;;; We need to play a trick to keep track of the number of nodes in the square
;;; between calls to MAKE-VERTEX-NODES and TESSELATE in MAKE-ELEMENT; we should
;;; find some way to restructure MAKE-ELEMENT so that this is not necessary.

(define square-domain-constructor
  (let ((width 0)
	(height 0)
	(h 0.)
	(k 0.))

    (list
     (lambda (m n)
       (set! width (+ m 2))
       (set! height (+ n 2))
       (let ((nodes (make-vector (* width height))))
	 (set! h (exact->inexact (/ (+ m 1))))
	 (set! k (exact->inexact (/ (+ n 1))))
	 (write-line `(h/k = ,(/ h k)))

	 (do ((i 0 (+ i 1)))
	     ((>= i width) nodes)
	   (do ((j 0 (+ j 1)))
	       ((>= j height))
	     (vector-set! nodes (+ (* i height) j)
			  (make-node (* i h) (* j k)))))))

     (lambda (nodes)

       ;; Triangulate:

       (let ((height-1 (- height 1))
	     (width-1 (- width 1))
	     (get-node (lambda (i) (vector-ref nodes i))))
	 (let column ((i 0) (results '(()())))
	   (if (< i width)
	       (column
		(+ i 1)
		(let row ((j 0) (edges (car results)) (faces (cadr results)))
		  (if (< j height)
		      (let ((me (+ (* i height) j))
			    (north (+ (* i height) j 1))
			    (east (+ (* i height) j height))
			    (ne (+ (* i height) j height 1)))

			(if (= j height-1)
			    (list edges faces)
			    (if (= i width-1)
				(row (+ j 1)
				     (cons (map get-node (list me north))
					   edges)
				     faces)
				(let ((f1 (list	(vector-ref nodes me)
						(vector-ref nodes north)
						(vector-ref nodes ne)))
				      (f2 (list (vector-ref nodes me)
						(vector-ref nodes ne)
						(vector-ref nodes east))))
				  (row (+ j 1)
				       (append (map
						(lambda (l)
						  (map get-node l))
						(list (list me north)
						      (list north ne)
						      (list me ne)
						      (list ne east)
						      (list me east)))
					       edges)
				       (cons f1 (cons f2 faces))))))))))
	       results))))

     (lambda () (list h k)))))

(define make-square-domain-vertices
  (car square-domain-constructor))

(define square-domain-triangulation
  (cadr square-domain-constructor))

(define square-domain-element-size
  (caddr square-domain-constructor))

(define (dirichlet-boundary? node)
  (let ((x (node:get-x node))
	(y (node:get-y node)))
    (memq #t (map almost-zero? (list x y (- 1 x) (- 1 y))))))

(define (cauchy-boundary? node)
  (let ((x (node:get-x node))
	(t (node:get-y node))
	(threshold (* .75 (cadr (square-domain-element-size)))))
    (or (memq #t (map almost-zero? (list x (- 1 x))))
	(< t threshold))))


;;; Square domain with right triangle attached:

(define (hat-boundary? node)
  (let* ((x (node:get-x node))
	 (t (node:get-y node))
	 (sizes (square-domain-element-size))
	 (dx (car sizes))
	 (dt (cadr sizes))
	 (m (/ dt dx))
	 (b (- 1 m)))

    (or (< t (* .75 dt))
	(> x (- 1 (* .25 dx)))
	(> t (- (+ (* x m) b) (* .25 dt)))
	(< x (* .25 dx)))))


;;; Random domains:

(define (make-random-square-domain-vertices n)
  (let* ((border (inexact->exact (floor (sqrt n))))
         (size (+ n (* 4 border) 4))
         (nodes (make-vector size)))

    (write-line `(border: ,(+ (* 4 border) 4) total: ,size))

    ;; Create the border (top and bottom):

    (let ((h (/ 1. (+ border 1))))
      (let ((border+2 (+ border 2)))
        (do ((i 0 (+ i 1)))
            ((>= i border+2))
          (let ((new-node (make-node (* i h) 0.)))
            (vector-set! nodes (* 2 i) new-node))
          (let ((new-node (make-node (* i h) 1.)))
            (vector-set! nodes (+ (* 2 i) 1) new-node))))

      ;; Left and right:

      (let ((2*border+2 (+ (* 2 border) 2)))
        (do ((j 1 (+ j 1)))
            ((> j border))
          (let ((new-node (make-node 0. (* j h))))
            (vector-set! nodes (+ (* 2 j) 2*border+2) new-node))
          (let ((new-node (make-node 1. (* j h))))
            (vector-set! nodes (+ (* 2 j) 2*border+2 1) new-node))))

      ;; Make internal nodes:

      (do ((i (+ (* 4 border) 4) (+ i 1)))
          ((>= i size) nodes)
        (let* ((x (random 1.))
               (y (random 1.))
               (new-node (make-node x y)))
          (vector-set! nodes i new-node))))))

(define (make-random-domain-vertices n)
  (let ((nodes (make-vector n)))
    (do ((i 0 (+ i 1)))
	((>= i n) nodes)
      (let* ((x (random 1.))
	     (y (random 1.))
	     (new-node (make-node x y #f)))
	(vector-set! nodes i new-node)))))

(define (random-domain-triangulation nodes)

  ;; Triangulate:

  (let ((n (vector-length nodes)))
    (do ((i 0 (+ i 1)))
	((>= i n))
      (node:set-boundary! (vector-ref nodes i) #f)))

  (let ((chull (convex-hull nodes)))
    (write-line `(the convex hull has ,(length chull) nodes...))
    (for-each
     (lambda (e) (node:set-boundary! (org e) #t))
     chull))

  (list (map (lambda (e)
	       (list (org e) (dest e)))
	     (list-edges))
	(map (lambda (f)
	       (map org f))
	     (list-faces))))

(define (make-not-so-random-domain-vertices n max-r)
  (let ((nodes (make-vector n)))
    (do ((i 0 (+ i 1)))
        ((>= i n) nodes)
      (let try ((x (random 1.)) (y (random 1.)))
        (let ((new-node (make-node x y #f)))
          (let loop ((j 0))
            (if (< j i)
                (if (< (nodal-distance new-node (vector-ref nodes j)) max-r)
                    (try (random 1.) (random 1.))
                    (loop (+ j 1)))
                (vector-set! nodes i new-node))))))))


;;; Making triangular domains:

(define triangular-domain-constructor
  (let ((dx 0.)
	(dy 0.))
    (list

     (lambda (n)
       ;; N is the number of nodes along the base.
       (let ((nodes (make-vector (/ (* n (+ n 1)) 2))))

	 (set! dx (/ 1. (- n 1)))
	 (set! dy dx)

	 (write-line `(making ,(* n (- (* 2 n) 1)) nodes...))

	 (let ((count 0))
	   (do ((i 0 (+ i 1))
		(y 0. (+ y dy))
		(row (- n 1) (- row 1)))
	       ((>= i n) nodes)
	     (do ((j 0 (+ j 1))
		  (x (/ y 2) (+ x dx)))
		 ((> j row))
	       (vector-set! nodes count (make-node x y))
	       (set! count (+ count 1)))))))

     (lambda () (list dx dy)))))

(define make-triangular-domain-vertices (car triangular-domain-constructor))
(define triangular-domain-element-size (cadr triangular-domain-constructor))

(define (triangular-boundary? node)
  (let ((dt (cadr (triangular-domain-element-size)))
	(x (node:get-x node))
	(t (node:get-y node)))
    (or (memq #t (map almost-zero? (list (- t (* 2 x)) (- t (- 2 (* 2 x))))))
	(< t (* .75 dt)))))


;;; Construct a true hat domain:

(define *hat-vertices-data* '())

(define (make-hat-vertices t-count x-count ratio)

  ;; Build a "hat domain," where the rectangular part is the unit square and
  ;; contains X-COUNT by T-COUNT nodes.  RATIO is the slope of the leg of the
  ;; hat.

  (let* ((node-count (+ (* t-count x-count) (/ (* x-count (- x-count 1)) 2)))
	 (nodes (make-vector node-count))
	 (dx (/ 1. (- x-count 1)))
	 (dt (/ 1. (- t-count 1))))

    (write-line `(,node-count nodes))

    (set! *hat-vertices-data* (list dx dt ratio))

    ;; The square part:

    (do ((i 0 (+ i 1)))
	((>= i x-count))
      (do ((j 0 (+ j 1)))
	  ((>= j t-count))
	(let ((x (* i dx))
	      (t (* j dt)))
	  (vector-set! nodes (+ (* i t-count) j) (make-node x t)))))

    ;; The triangle:

    (let ((count (* t-count x-count)))
      (set! dt (* ratio dx))

      (do ((j 1 (+ j 1)))
	  ((>= j x-count))
	(let* ((t (+ 1 (* j dt)))
	       (x0 (/ (- t 1) ratio 2)))
	  (do ((i 0 (+ i 1)))
	      ((>= i (- x-count j)))
	    (let ((x (+ x0 (* i dx))))
	      (vector-set! nodes count (make-node x t))
	      (set! count (+ count 1)))))))
    nodes))

(define (true-hat-boundary? node)
  (let ((dx (car *hat-vertices-data*))
	(dt (cadr *hat-vertices-data*))
	(ratio (caddr *hat-vertices-data*)))
    (let ((x (node:get-x node))
	  (t (node:get-y node))
	  (dx/4 (/ dx 4))
	  (3/4*dt (* 3/4 dt))
	  (dx/8 (/ dx 8)))
      (or (and (<= t (+ 1 dx/4))
	       (or (<= t 3/4*dt)
		   (<= x dx/4)
		   (>= x (- 1 dx/4))))
	  (and (>= t 1)
	       (or (and (<= x (+ .5 dx/8))
			(<= x (+ (/ (- t 1) ratio 2) dx/8)))
		   (and (>= x (- .5 dx/8))
			(>= x (- (- 1 (/ (- t 1) ratio 2)) dx/8)))))))))
