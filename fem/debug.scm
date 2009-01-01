;;; For debugging purposes, this program draws the output of the Delaunay
;;; triangulation program, TRIANGULATE, in delaunay.scm.
;;;
;;; We assume that the sites are subsets of the unit square [0,1]x[0,1].

(declare (usual-integrations))


;;; Draw the mesh:

(define *delaunay-device* 'undefined)
(define *debugging-info* '())

(define draw
  (let ((background "black")
	(cursor "white")
	(line "blue")
	(boundary "red")
	(node "white")
	(boundary-node "purple"))
    (lambda (nodes)
      (if (eq? *delaunay-device* 'undefined)
	  (set! *delaunay-device* (make-graphics-device 'x))
	  (graphics-clear *delaunay-device*))

      (let* ((dev *delaunay-device*)
	     (b .02)
	     (elist (set-coordinate-limits dev nodes b .005))
	     (ex (car elist))
	     (ey (cadr elist)))
	(graphics-operation dev 'set-foreground-color line)
	(graphics-operation dev 'set-background-color background)
	(graphics-operation dev 'set-mouse-color cursor)
	(graphics-clear dev)

	(for-each
	 (lambda (e)
	   (let* ((org-e (car e))
		  (dest-e (cadr e))
		  (org (node:get-coords org-e))
		  (dest (node:get-coords dest-e))
		  (org-boundary? (node:boundary? org-e))
		  (dest-boundary? (node:boundary? dest-e)))

	     (if (and org-boundary? dest-boundary?)
		 (graphics-operation dev 'set-foreground-color boundary))

	     (graphics-move-cursor dev (vector-ref org 0) (vector-ref org 1))
	     (graphics-drag-cursor dev (vector-ref dest 0) (vector-ref dest 1))


	     (if (and org-boundary? dest-boundary?)
		 (graphics-operation dev 'set-foreground-color line))))

	 *debugging-info*)

	(graphics-operation dev 'set-foreground-color node)

	(for-each
	 (lambda (n)
	   (let ((x (node:get-x n))
		 (y (node:get-y n)))

	     (if (node:boundary? n)
		 (graphics-operation dev 'set-foreground-color boundary-node))

	     (graphics-draw-line dev (- x ex) (- y ey) (+ x ex) (+ y ey))
	     (graphics-draw-line dev (- x ex) (+ y ey) (+ x ex) (- y ey))

	     (if (node:boundary? n)
		 (graphics-operation dev 'set-foreground-color node))))

	 (vector->list nodes))))))

(define (set-coordinate-limits dev nodes border edge)
  (apply
   (lambda (x-left y-bottom x-right y-top)

     (if (= x-left x-right)
	 (begin
	   (set! x-left (- x-left .5))
	   (set! x-right (+ x-right .5))))

     (if (= y-top y-bottom)
	 (begin
	   (set! y-bottom (- y-bottom .5))
	   (set! y-top (+ y-top .5))))

     (let ((dx (- x-right x-left))
	   (dy (- y-top y-bottom)))
       (graphics-set-coordinate-limits
	dev (- x-left (* dx border)) (- y-bottom (* dy border))
	(+ x-right (* dx border)) (+ y-top (* dy border)))
       (list (* edge dx) (* edge dy))))
   (bounding-box (vector->list nodes) node:get-coords)))

(define (close)
  (if (not (eq? *delaunay-device* 'undefined))
      (begin
	(graphics-close *delaunay-device*)
	(set! *delaunay-device* 'undefined))))


;;; Quick and easy way to dump a vector, matrix, whatever into a file:

(define (dump obj file-name)
  ;; Compound objects should always come before primitive ones, because they
  ;; may be implemented from primitives like LISTs or VECTORs.
  (cond ((matrix? obj)
         (write-line `(dumping matrix to file ,file-name))
         (let ((port (open-output-file file-name)))
           (print-matrix obj port)
           (close-output-port port)))

        ((sparse-matrix? obj)
         (write-line `(dumping sparse-matrix to file ,file-name))
         (let ((port (open-output-file file-name)))
           (print-sparse-matrix obj port)
           (close-output-port port)))

        ((vector? obj)
         (write-line `(dumping vector to file ,file-name))

         (let ((n (vector-length obj))
               (port (open-output-file file-name)))

           (do ((i 0 (+ i 1)))
               ((>= i n))
             (display (vector-ref obj i) port)
             (newline port))

           (close-output-port port)))

        ((list? obj)
         (write-line `(dumping list to file ,file-name))
         (let ((port (open-output-file file-name)))

           (for-each
            (lambda (x)
              (display x port)
              (newline port))
            obj)

           (close-output-port port)))

        (else
         (error "Object must be a VECTOR, MATRIX, or LIST -- DUMP"))))


;;; Dump a matrix into a Maple-readable file:

(define matrix->maple
  (let ((variable-name "foo"))
    (lambda (A file-name)
      (let* ((port (open-output-file file-name))
	     (print (lambda (obj) (display obj port))))

	(print (string-append variable-name " := ["))
	(print #\newline)

	(let ((m (matrix-row-count A))
	      (n (matrix-column-count A)))
	  (print (string-append "[" (number->string (matrix-ref A 0 0))))

	  (do ((j 1 (+ j 1)))
	      ((>= j n))
	    (print (string-append "," (number->string (matrix-ref A 0 j)))))

	  (print #\])

	  (do ((i 1 (+ i 1)))
	      ((>= i m))
	    (print (string-append (list->string '(#\, #\newline #\[))
				  (number->string (matrix-ref A i 0))))

	    (do ((j 1 (+ j 1)))
		((>= j n))
	      (print (string-append "," (number->string (matrix-ref A i j)))))

	    (print #\]))

	  (print "]:")
	  (print #\newline)
	  (print (string-append variable-name
				" := array("
				variable-name
				"):"))
	  (print #\newline)
	  (print (string-append variable-name
				" := [seq([seq([i, j, "
				variable-name
				"[i, j]], i=1.."
				(number->string m)
				") ], j=1.."
				(number->string n)
				")]:"))
	  (print #\newline)
	  (print (string-append "plots[surfdata]("
				variable-name
				", axes=frame, style=wireframe);"))
	  (print #\newline))

	(close-output-port port)))))


;;; Compute the error vector:

(define (compute-error nodes v f)
  ;; V should be the output of SOR.
  ;; F computes the initial condition/solution, given a node.
  (let ((e (make-matrix (vector-length v) 3))
        (size (vector-length nodes))
        (max 0)
        (max-index 0))
    (let loop ((i 0) (j 0))
      (if (< j size)
          (if (node:boundary? (vector-ref nodes j))
              (loop i (+ j 1))
              (begin
                (matrix-set! e i 0 (f (vector-ref nodes j)))
                (matrix-set! e i 1 (vector-ref v i))
                (matrix-set! e i 2 (abs (- (f (vector-ref nodes j))
                                           (vector-ref v i))))
                (if (>= (matrix-ref e i 2) max)
                    (begin
                      (set! max (matrix-ref e i 2))
                      (set! max-index i)))
                (loop (+ i 1) (+ j 1))))

          (begin
            (write-line `(maximum error: ,max at node ,max-index))
            e)))))


;;; Generate results that we can plot with MATLAB:

(define (compute-results nodes v f)

  ;; V should be the output of SOR.

  (let ((e (make-matrix (interior-node-count nodes) 4))
        (size (vector-length nodes)))
    (let loop ((i 0) (j 0))
      (if (< j size)
          (if (node:boundary? (vector-ref nodes j))
	      (loop i (+ j 1))
	      (let ((node (vector-ref nodes j)))
		(matrix-set! e i 0 (node:get-x node))
		(matrix-set! e i 1 (node:get-y node))
		(matrix-set! e i 2 (vector-ref v i))
		(matrix-set! e i 3 (f node))
		(loop (+ i 1) (+ j 1))))
          e))))

(define (interior-node-count nodes)
  (let ((n (vector-length nodes)))
    (let loop ((count 0) (i 0))
      (if (< i n)
	  (if (node:boundary? (vector-ref nodes i))
	      (loop count (+ i 1))
	      (loop (+ count 1) (+ i 1)))
	  count))))

;;; Save the data back to the nodes:

(define (store-results! nodes v)
  ;; V should be the output of SOR.
  (let ((size (vector-length nodes)))
    (let loop ((i 0) (j 0))
      (if (< j size)
	  (if (node:boundary? (vector-ref nodes j))
	      (loop i (+ j 1))
	      (begin
		(node:set-value! (vector-ref nodes j) (vector-ref v i))
		(loop (+ i 1) (+ j 1))))
	  (write-line `(,i interior nodes))))))
