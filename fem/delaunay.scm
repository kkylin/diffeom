;;; This implements the Delaunay triangulation algorithm described in:
;;;
;;; "Primitives for the manipulation of general subdivisions and the
;;;  computation of Voronoi diagrams," ACM transactions on graphics, Vol. 4,
;;; No. 2, April 1985, P. 74-123.
;;; Leonidas Guibas and Jorge Stolfi, Xerox PARC and Stanford University.

;;; There is something better now; take a look at the Geometry Center's home
;;; page at "http://www.geom.umn.edu/", and check out the quickhull algorithm.

(declare (usual-integrations))


;;; This is kind of useful for the FEM stuff:

(define (delaunay-triangulation nodes)
  (triangulate nodes)
  (list (map (lambda (e)
	       (list (org e) (dest e)))
	     (list-edges))
	(map (lambda (f)
	       (map org f))
	     (list-faces))))


;;; This is the divide-and-conquer algorithm.  The nodal data structure should
;;; provide the methods NODE:GET-X and NODE:GET-Y.

(define (triangulate nodes)
  ;; NODES should be a vector containing nodal data structures.
  (set! *delaunay-edges* (make-dynamic-table))
  (delaunay (sort (vector->list nodes) lexicographic<)))

(define (lexicographic< n1 n2)
  (let ((delta-x (- (node:get-x n2) (node:get-x n1))))
    (if (almost-zero? delta-x)
	(< (node:get-y n1) (node:get-y n2))
	(> delta-x 0))))

(define (delaunay S)
  ;; This returns the counterclockwise convex hull edge out of
  ;; the leftmost vertex and the clockwise convex hull edge out
  ;; of the rightmost vertex.
  ;;
  ;; S is assumed to be a list of nodes, sorted along the abscissa.

  (cond
   ((< (length S) 2)
    (error "Need at least two nodes to triangulate -- DELAUNAY"))
   ((= (length S) 2)
    ;; Let s1, s2 be the two sites, in sorted order, and create an edge from s1
    ;; to s2:
    (let* ((s1 (car S))
           (s2 (cadr S))
           (a (make-edge)))
      (set-org! a s1)
      (set-dest! a s2)
      (list a (sym a))))

   ((= (length S) 3)
    ;; Let s1, s2, and s3 be the three sites, in sorted order.
    ;; Create edges a connecting s1 to s2 and b connecting s2 to s3:
    (let* ((s1 (car S))
           (s2 (cadr S))
           (s3 (caddr S))
           (a (make-edge))
           (b (make-edge)))

      (splice (sym a) b)
      (set-org! a s1)
      (set-dest! a s2)
      (set-org! b s2)
      (set-dest! b s3)

      ;; Now close the triangle:
      (cond
       ((ccw s1 s2 s3)
        (connect b a)
        (list a (sym b)))

       ((ccw s1 s3 s2)
        (let ((c (connect b a)))
          (list (sym c) c)))

       (else ; the three points are collinear
        (list a (sym b))))))

   (else
    ;; |S| > 3.  Let L and R be the left and right halves of S.
    (let* ((L&R (halve S))
           (L (car L&R))
           (R (cadr L&R))

           (ldo&ldi (delaunay L))
           (ldo (car ldo&ldi))
           (ldi (cadr ldo&ldi))

           (rdi&rdo (delaunay R))
           (rdi (car rdi&rdo))
           (rdo (cadr rdi&rdo)))

      ;; Compute the lower common tangent of L and R:

      (call-with-current-continuation
        (lambda (exit)
          (do () (#f)
            (cond ((left-of (org rdi) ldi)
                   (set! ldi (lnext ldi)))
                  ((right-of (org ldi) rdi)
                   (set! rdi (rprev rdi)))
                  (else (exit 'done))))))

      ;; Create a first cross edge basel from rdi.Org to ldi.Org:

      (let* ((basel (connect (sym rdi) ldi))
             (valid (lambda (e) (right-of (dest e) basel))))

        (if (node= (org ldi) (org ldo)) (set! ldo (sym basel)))
        (if (node= (org rdi) (org rdo)) (set! rdo basel))

        ;; This is the merge loop:
        ;; Locate the first L point (lcand.Dest) to be encountered by the
        ;; rising bubble, and delete L edges out of basel.Dest that fail the
        ;; circle test.

        (call-with-current-continuation
         (lambda (exit)
           (do () (#f)

             (let ((lcand (onext (sym basel))))
               (if (valid lcand)
                   (do ()
                       ((not (in-circle (dest basel) (org basel) (dest lcand)
                                        (dest (onext lcand)))))
                     (let ((t (onext lcand)))
                       (delete-edge lcand)
                       (set! lcand t))))

               ;; Symmetrically, locate the first R point to be hit, and delete
               ;; R edges:

               (let ((rcand (oprev basel)))
                 (if (valid rcand)
                     (do ()
                         ((not (in-circle (dest basel) (org basel) (dest rcand)
                                          (dest (oprev rcand)))))
                       (let ((t (oprev rcand)))
                         (delete-edge rcand)
                         (set! rcand t))))

                 ;; If both lcand and rcand are invalid, then basel is the
                 ;; upper common tangent:

                 (if (and (not (valid lcand)) (not (valid rcand)))
                     (exit 'done))

                 ;; The next cross edge is to be connected to either lcand.Dest
                 ;; or rcand.Dest.  If both are valid, then choose the
                 ;; appropriate one using the InCircle test:

                 (if (or (not (valid lcand))
                         (and (valid rcand)
                              (in-circle (dest lcand) (org lcand) (org rcand)
                                         (dest rcand))))
                     ;; Add cross edge basel from rcand.Dest to basel.Dest:
                     (set! basel (connect rcand (sym basel)))
                     ;; Else add cross edge basel from basel.Org to lcand.Dest:
                     (set! basel (connect (sym basel) (sym lcand))))))))))

      (list ldo rdo)))))


;;; Miscellaneous functions:

(define (halve l)
  ;; Split L down the middle
  (let loop ((hcount 0) (tcount (length l)) (head '()) (tail l))
    (if (<= tcount hcount)
        (list (reverse head) tail)
        (loop (+ hcount 1) (- tcount 1) (cons (car tail) head) (cdr tail)))))

(define (node= n1 n2)
  (apply equal? (map node:get-coords `(,n1 ,n2))))


;;; Use the above algorithm to compute the convex hull:

(define (convex-hull nodes)
  (let* ((e (cadr (triangulate nodes)))
	 (e.org (org e)))
    (let loop ((l (list e)))
      (let ((next (lnext (car l))))
	(if (node= (org next) e.org)
	    l
	    (loop (cons next l)))))))
