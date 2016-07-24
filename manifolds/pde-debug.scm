#|

This file is part of DIFFEOM, a system for solving
differential equations on manifolds.

Copyright (C) 2016 by Kevin K Lin
<kkylin@alum.mit.edu>

This program is free software; you can redistribute
it and/or modify it under the terms of the GNU
General Public License as published by the Free
Software Foundation; either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General
Public License along with this program; if not, write
to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301 USA.

|#

;;; This file is based on fem/debug.scm.

(declare (usual-integrations))


;;; Drawing a domain requires some care with element edges:

(define (draw-domain domain)
  (draw-charts-in-domain (manifold:get-finite-atlas domain) domain))


;;; Drawing a chart is much simpler:

(define (draw-chart chart)
  (let ((edges (complex->edges (chart:get-complex chart)))
	(nodes (chart:get-nodes chart)))
    (apply 2d-draw
	   (append
	    (list (chart:get-nodes chart) edges node:get-coords draw-line)
	    (bounding-box nodes node:get-coords)))))

(define (draw-chart-in-domain chart domain)
  (draw-charts-in-domain (list chart) domain))

(define (draw-charts-in-domain charts domain)
  (let ((edges (append-map (compose complex->edges chart:get-complex) charts))
	(nodes (append-map chart:get-nodes charts)))
    (apply 2d-draw
	   (append
	    (list nodes edges node:get-point draw-edge-in-chart)
	    (bounding-box (manifold:get-nodes domain) node:get-point)))))


;;; Actually draw something:

(define *planar-device* 'undefined)

(define (make-2d-draw colors)
  (let ((background (vector-ref colors 0))
	(cursor (vector-ref colors 1))
	(line (vector-ref colors 2))
	(boundary (vector-ref colors 3))
	(loc-bound (vector-ref colors 4))
	(node (vector-ref colors 5))
	(bnode (vector-ref colors 6))
	(gnode (vector-ref colors 7))

	(border-frac .05)
	(cross-frac .005))

    (lambda (nodes edges get-coords draw-edge
		   x-left y-bottom x-right y-top)

      ;; Report the dimensions:

      (write-line `(x range: ,x-left to ,x-right))
      (write-line `(y range: ,y-bottom to ,y-top))

      ;; Open a graphics window, if there isn't one already:

      (if (eq? *planar-device* 'undefined)
	  (set! *planar-device* (make-graphics-device 'x))
	  (graphics-clear *planar-device*))

      (let ((dev *planar-device*)
	    (bx (* border-frac (- x-right x-left)))
	    (by (* border-frac (- y-top y-bottom)))
	    (ex (* cross-frac (- x-right x-left)))
	    (ey (* cross-frac (- y-top y-bottom))))

	(if (zero? bx) (set! bx .5))
	(if (zero? by) (set! by .5))
	(if (zero? ex) (set! ex cross-frac))
	(if (zero? ey) (set! ey cross-frac))

	;; Set up the window:

	(write-line '(setting up...))

	(graphics-set-coordinate-limits
	 dev (- x-left bx) (- y-bottom by) (+ x-right bx) (+ y-top by))
	(graphics-operation dev 'set-foreground-color line)
	(graphics-operation dev 'set-background-color background)
	(graphics-operation dev 'set-mouse-color cursor)
	(graphics-clear dev)
	(graphics-enable-buffering dev)

	;; First, draw the edges (the graphics are cooler if one sorts the
	;; edges first):

	(write-line `(drawing ,(length edges) edges...))

	(let ((color 'line))
	  (for-each
	   (lambda (e)
	     (let* ((org-e (car e))
		    (dest-e (cadr e))
		    (org (get-coords org-e))
		    (dest (get-coords dest-e)))

	       (cond ((and (node:boundary? org-e) (node:boundary? dest-e))

		      (if (not (eq? color 'boundary))
			  (begin
			    (set! color 'boundary)
			    (graphics-operation
			     dev 'set-foreground-color boundary))))

		     ((and (node:local-boundary? org-e)
			   (node:local-boundary? dest-e))

		      (if (not (eq? color 'local-boundary))
			  (begin
			    (set! color 'local-boundary)
			    (graphics-operation
			     dev 'set-foreground-color loc-bound))))

		     (else

		      (if (not (eq? color 'line))
			  (begin
			    (set! color 'line)
			    (graphics-operation
			     dev 'set-foreground-color line)))))

	       (draw-edge dev e)))

	   edges))

	;; Next, draw the nodes:

	(write-line `(drawing ,(length nodes) nodes...))

	(let ((color 'node))
	  (graphics-operation dev 'set-foreground-color node)

	  (for-each
	   (lambda (n)
	     (let* ((coords (get-coords n))
		    (x (vector-ref coords 0))
		    (y (vector-ref coords 1)))

	       (cond ((node:boundary? n)

		      (if (not (eq? color 'bnode))
			  (begin
			    (set! color 'bnode)
			    (graphics-operation
			     dev 'set-foreground-color bnode)))

		      (draw-boundary-node dev x y ex ey))

		     ((node:constrained? n)

		      (if (not (eq? color 'gnode))
			  (begin
			    (set! color 'gnode)
			    (graphics-operation
			     dev 'set-foreground-color gnode)))

		      (draw-glued-node dev x y ex ey))

		     (else

		      (if (not (eq? color 'node))
			  (begin
			    (set! color 'node)
			    (graphics-operation
			     dev 'set-foreground-color node)))

		      (draw-node dev x y ex ey)))))
	   nodes))

	(write-line '(flushing buffers...))
	(graphics-disable-buffering dev))
      'done)))


;;; Different color schemes for different purposes:

(define *standard-colors*
  (vector "black"          ;; Background.
	  "white"          ;; Cursor.
	  "blue"           ;; Regular edges.
	  "red"            ;; Boundary edges.
	  "purple"         ;; Local boundaries.
	  "white"          ;; Regular node color.
	  "orange"         ;; Boundary node color.
	  "green"))        ;; "Glued" node color.

(define *boring-colors*
  (vector "white"          ;; Background.
	  "black"          ;; Cursor.
	  "black"          ;; Regular edges.
	  "black"          ;; Boundary edges.
	  "black"          ;; Local boundaries.
	  "black"          ;; Regular node color.
	  "black"          ;; Boundary node color.
	  "black"))        ;; "Glued" node color.

(define *print-colors*
  (vector "white"          ;; Background.
	  "black"          ;; Cursor.
	  "blue"           ;; Regular edges.
	  "red"            ;; Boundary edges.
	  "yellow"         ;; Local boundaries.
	  "black"          ;; Regular node color.
	  "red"            ;; Boundary node color.
	  "black"))        ;; "Glued" node color.

(define 2d-draw (make-2d-draw *standard-colors*))


;;; Methods for drawing nodes:

(define (draw-node dev x y ex ey)

  ;; An 'x' is a regular node:

  (graphics-draw-line dev (- x ex) (- y ey) (+ x ex) (+ y ey))
  (graphics-draw-line dev (- x ex) (+ y ey) (+ x ex) (- y ey)))

(define (draw-boundary-node dev x y ex ey)

  ;; A square is a boundary node:

  (let ((x-ex (- x ex))
	(x+ex (+ x ex))
	(y-ey (- y ey))
	(y+ey (+ y ey)))

    (graphics-move-cursor dev x-ex y-ey)
    (graphics-drag-cursor dev x-ex y+ey)
    (graphics-drag-cursor dev x+ex y+ey)
    (graphics-drag-cursor dev x+ex y-ey)
    (graphics-drag-cursor dev x-ex y-ey)))

(define (draw-glued-node dev x y ex ey)

  ;; A triangle is a node glued to another chart:

  (let ((x-ex (- x ex))
	(x+ex (+ x ex))
	(y-ey (- y ey))
	(y+ey (+ y ey)))

    (graphics-draw-line dev x-ex y-ey x y+ey)
    (graphics-draw-line dev x y+ey x+ex y-ey)
    (graphics-draw-line dev x+ex y-ey x-ex y-ey)))


;;; Getting rid of the window is sometimes useful, too!

(define (close)
  (if (not (eq? *planar-device* 'undefined))
      (begin
	(graphics-close *planar-device*)
	(set! *planar-device* 'undefined))))


;;; In order to get edges drawn correctly, we need to walk around the
;;; parametrized path in the chart and then map it into the manifold (the
;;; straight line in the manifold seldom works correctly):

(define draw-edge-in-chart
  (let* ((step-count 40)
	 (dt (exact->inexact (/ step-count))))
    (lambda (dev edge)
      (let* ((org (node:get-coords (car edge)))
	     (dir (vector:- (node:get-coords (cadr edge)) org))
	     (f (chart:get-inverse-map (node:get-chart (car edge)))))
	(let loop ((i 0) (p (f org)) (q (f (vector:+ (vector:* dt dir) org))))
	  (if (< i step-count)
	      (let ((i+1 (+ i 1)))
		(graphics-draw-line dev
				    (vector-ref p 0) (vector-ref p 1)
				    (vector-ref q 0) (vector-ref q 1))
		(loop i+1 q (f (vector:+ (vector:* (* i+1 dt) dir) org))))
	      'done))))))

;;; Draw a straight line:

(define (draw-line dev edge)
  (let ((org (node:get-coords (car edge)))
	(dest (node:get-coords (cadr edge))))
    (graphics-draw-line dev
			(vector-ref org 0) (vector-ref org 1)
			(vector-ref dest 0) (vector-ref dest 1))))


;;; What we have so far actually doesn't work!  Even though we get more
;;; equations than unknowns, the rank is smaller than the number of unknowns
;;; (i.e. the system is underdetermined).  What can be going wrong?

;;; 1. Something is wrong in applying FEM, most notably FEM-DISCRETIZE in
;;;    pde-tools.scm.

;;; 2. Something is wrong in assembling the equations after that, as in
;;;    COMBINE-EQUATIONS.

;;; 3. Something is wrong in generating constraints.

;;; 4. Something is wrong in collecting the equations and constraints.

;;; 5. Something is wrong in *controlling* when to generate constraints.

;;; (5) is the hardest one to fix, but it's most likely what's wrong with the
;;; approach.

;;; We can do a end-to-end test of 1, 2, and 4 by extracting the local matrices
;;; and comparing them with the local results.  (There shouldn't be any bugs in
;;; FEM, right?)  It turns out that the specific example we have does pass the
;;; test.  So we need to figure out if there are problems with the constraints
;;; (items 3 and 5).


(define (extract-local-matrix mat chart)
  (let ((icount 0)
	(ncount-1 (- (matrix-column-count mat) 1))
	(result #f)
	(indices '()))

    (let loop ((nodes (chart:get-nodes chart)) (count 0) (ilist '()))
      (if (null? nodes)
	  (begin
	    (set! icount count)
	    (set! result (make-matrix icount (+ icount 1)))
	    (set! indices (sort ilist <)))
	  (let ((id (node:get-id (car nodes))))
	    (if (number? id)
		(loop (cdr nodes) (+ count 1) (cons id ilist))
		(loop (cdr nodes) count ilist)))))

    (do ((ilist indices (cdr ilist))
	 (m 0 (+ m 1)))
	((null? ilist) result)
      (let ((i (car ilist)))
	(do ((jlist indices (cdr jlist))
	     (n 0 (+ n 1)))
	    ((null? jlist))
	  (let ((j (car jlist)))
	    (matrix-set! result m n (matrix-ref mat i j))))
	(matrix-set! result m icount (matrix-ref mat i ncount-1))))))

(define (compute-local-matrix source chart)
  (sparse->matrix (assemble-equations
		   source (list->vector (chart:get-nodes chart)))))

(define (compare-matrices A B)
  (let ((m (matrix-row-count A))
	(n (matrix-column-count A)))
    (do ((i 0 (+ i 1)))
	((>= i m))
      (do ((j 0 (+ j 1)))
	  ((>= j n))
	(let ((diff (- (matrix-ref A i j) (matrix-ref B i j))))
	  (if (not (zero? diff))
	      (write-line
	       `((,i ,j)
		 (- ,(matrix-ref A i j) ,(matrix-ref B i j))
		 => ,diff))))))))


;;; Test if a column of a matrix is all almost zero.

(define (matrix:null-column? mat j)
  (let ((m (matrix-row-count mat)))
    (let loop ((i 0))
      (if (< i m)
	  (if (almost-zero? (matrix-ref mat i j))
	      (loop (+ i 1))
	      #f)
	  #t))))

(define (matrix:null-row? mat i)
  (let ((n (matrix-column-count mat)))
    (let loop ((j 0))
      (if (< j n)
	  (if (almost-zero? (matrix-ref mat i j))
	      (loop (+ j 1))
	      #f)
	  #t))))

(define (matrix:find-null-columns mat)
  (let ((n (matrix-column-count mat)))
    (let loop ((j 0) (results '()))
      (if (< j n)
	  (if (matrix:null-column? mat j)
	      (loop (+ j 1) (cons j results))
	      (loop (+ j 1) results))
	  results))))

(define (matrix:find-null-rows mat)
  (let ((m (matrix-row-count mat)))
    (let loop ((i 0) (results '()))
      (if (< i m)
	  (if (matrix:null-row? mat i)
	      (loop (+ i 1) (cons i results))
	      (loop (+ i 1) results))
	  results))))


;;; Compute some errors:

(define (process-absolute-errors process nodes f v)
  (apply process
	 (append-map
	  (lambda (node)
	    (let ((index (node:get-id node)))
	      (if (and (number? index) (not (node:constrained? node)))
		  (list (abs (- (f node) (vector-ref v index))))
		  '())))
	  nodes)))

(define (max-error domain f v)
  (process-absolute-errors max (manifold:get-nodes domain) f v))

(define (min-error domain f v)
  (process-absolute-errors min (manifold:get-nodes domain) f v))

(define (avg-error domain f v)
  (/ (process-absolute-errors + (manifold:get-nodes domain) f v)
     (vector-length v)))

(define (process-relative-errors process nodes f v)
  (apply process
	 (append-map
	  (lambda (node)
	    (let ((index (node:get-id node)))
	      (if (number? index)
		  (list (relative-error (vector-ref v index) (f node)))
		  '())))
	  nodes)))

(define (max-relative-error domain f v)
  (process-relative-errors max (manifold:get-nodes domain) f v))

(define (min-relative-error domain f v)
  (process-relative-errors min (manifold:get-nodes domain) f v))

(define (avg-relative-error domain f v)
  (/ (process-relative-errors + (manifold:get-nodes domain) f v)
     (vector-length v)))


;;; Somewhat more useful information about how bad the solutions are (if the
;;; number of nodes is reasonably small):

(define (display-results domain f v)
  (for-each
   (lambda (node)
     (let ((id (node:get-id node)))
       (if (number? id)
	   (write-line `((id = ,id)
			 (computed = ,(vector-ref v (node:get-id node)))
			 (actual = ,(f node))
			 (x = ,(node:get-x node) y = ,(node:get-y node)))))))
   (manifold:get-nodes domain)))



;;; This is so things can be drawn in MATLAB or Maple: Draw a rectangle that
;;; bounds the manifold, divide it into an MxN grid, drop the nodes and average
;;; the values.  Output is a matrix.

(define (manifold->grid m n domain f v diff)
  (let ((nodes (manifold:get-nodes domain)))

    (if (null? nodes)
	(error "Manifold contains no nodes! -- MANIFOLD->GRID"))

    ;; Save nodal values:

    (for-each
     (lambda (node)
       (let ((index (node:get-id node)))
	 (if (number? index)
	     (node:set-value! node (vector-ref v index)))))
     nodes)

    (let* ((p (node:get-point (car nodes)))
	   (x-min (vector-first p))
	   (x-max (vector-first p))
	   (y-min (vector-second p))
	   (y-max (vector-second p)))

      ;; Then find the bounding rectangle:

      (for-each
       (lambda (node)
	 (let* ((p (node:get-point node))
		(x (vector-first p))
		(y (vector-second p)))
	   (cond ((> x x-max) (set! x-max x))
		 ((< x x-min) (set! x-min x)))
	   (cond ((> y y-max) (set! y-max y))
		 ((< y y-min) (set! y-min y)))))
       (cdr nodes))

      ;; Next, create the matrix and start averaging:

      (let ((mat (make-matrix m n))
	    (count (make-matrix m n))
	    (dx (/ (- x-max x-min) m))
	    (dy (/ (- y-max y-min) n)))

	;; Collect the sums and count the number of nodes in each box:

	(for-each
	 (lambda (node)
	   (let* ((p (node:get-point node))
		  (x (- (vector-first p) x-min))
		  (y (- (vector-second p) y-min))
		  (i (inexact->exact (floor (/ x dx))))
		  (j (inexact->exact (floor (/ y dy)))))

	     ;; Silly fence-post conditions:

	     (if (= i m) (set! i (- m 1)))
	     (if (= j n) (set! j (- n 1)))

	     (matrix-set! count i j (+ 1 (matrix-ref count i j)))
	     (matrix-set! mat i j
			  (+ (diff (node:get-value node) (f node))
			     (matrix-ref mat i j)))))
	 nodes)

	;; Normalize:

	(do ((i 0 (+ i 1)))
	    ((>= i m) mat)
	  (do ((j 0 (+ j 1)))
	      ((>= j n))
	    (if (> (matrix-ref count i j) 0)
		(matrix-set!
		 mat i j
		 (/ (matrix-ref mat i j) (matrix-ref count i j))))))))))


;;; This lets save the results of a computation, if not the finite element
;;; stuff itself:

(define (node-states domain f v)
  (let ((nodes (manifold:get-nodes domain)))
    (if (null? nodes)
	#f
	(let* ((n (real-node-count nodes))
	       (mat (make-matrix n 4)))
	  (for-each
	   (lambda (node)
	     (let ((id (node:get-id node)))
	       (if (number? id)
		   (let ((p (node:get-point node)))
		     (matrix-set! mat id 0 (vector-first p))
		     (matrix-set! mat id 1 (vector-second p))
		     (matrix-set! mat id 2 (vector-ref v id))
		     (matrix-set! mat id 3 (f node))))))
	   nodes)
	  mat))))

(define (real-node-count nodes)
  (let loop ((count 0) (nodes nodes))
    (if (null? nodes)
	count
	(let ((id (node:get-id (car nodes))))
	  (if (number? id)
	      (loop (+ count 1) (cdr nodes))
	      (loop count (cdr nodes)))))))
