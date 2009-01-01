;;; This is a simple implementation of the quad-edge data structure and the
;;; associated edge functions on the edge algebra, as described in:
;;;
;;; "Primitives for the manipulation of general subdivisions and the
;;;  computation of Voronoi diagrams," ACM transactions on graphics, Vol. 4,
;;; No. 2, April 1985, P. 74-123.
;;; Leonidas Guibas and Jorge Stolfi, Xerox PARC and Stanford University.
;;;
;;; This is, in particular, meant to be used by the Delaunay triangulation
;;; algorithm.

(declare (usual-integrations))


;;; A global list of edges:

(define *delaunay-edges* (make-dynamic-table))


;;; The quad-edge data structure: This is the basic data structure
;;; representing the edges in an edge algebra.  It consists of a representative
;;; edge and its orbit under the operations Rot and Flip.

(define (make-edge-record id)
  (let ((orientation 0)
        (next (make-vector 4))
        ;; Auxiliary data fields, initialized to the null lists:
        (data (make-vector 4 'undefined))
        (mark (make-vector 4 'undefined)))

    ;; Main dispatcher:

    (define (me r)
      (lambda (msg)
        (case msg
          ((next) (vector-ref next r))
          ((set-next!) (lambda (val) (vector-set! next r val)))
          ((data) (vector-ref data r))
          ((set-data!) (lambda (val) (vector-set! data r val)))
          ((mark) (vector-ref mark r))
          ((set-mark!) (lambda (val) (vector-set! mark r val)))
          (else (error "Unknown request -- EDGE-RECORD")))))

    ;; Initialize the edges, which lie on a 2-sphere.

    (vector-set! next 0 (make-edge-ref me 0 orientation id))
    (vector-set! next 1 (make-edge-ref me 3 orientation id))
    (vector-set! next 2 (make-edge-ref me 2 orientation id))
    (vector-set! next 3 (make-edge-ref me 1 orientation id))

    me))


;;; An edge-reference is a triplet (e, r, f), where e is an edge record and r
;;; and f are the corresponding Rot and Flip degrees.

(define (make-edge-ref e r f id)
  (vector e r f id))

(define (get-edge-record e-ref)
  (vector-ref e-ref 0))

(define (get-rot-deg e-ref)
  (vector-ref e-ref 1))

(define (get-flip-deg e-ref)
  (vector-ref e-ref 2))

(define (get-edge-id e-ref)
  (vector-ref e-ref 3))

(define (get-edge-ring e)
  (let ((e.dest (dest e)))
    (let loop ((l (list e)) (e e))
      (let ((e.onext (onext e)))
        (if (node= e.dest (dest e.onext))
            l
            (loop (cons e.onext l) e.onext))))))


;;; Basic edge functions on which others are built:

(define (rot e-ref)
  (let ((e (get-edge-record e-ref))
        (r (get-rot-deg e-ref))
        (f (get-flip-deg e-ref))
        (id (get-edge-id e-ref)))
    (make-edge-ref e (modulo (+ r 1 (* 2 f)) 4) f id)))

(define (flip e-ref)
  (let ((e (get-edge-record e-ref))
        (r (get-rot-deg e-ref))
        (f (get-flip-deg e-ref))
        (id (get-edge-id e-ref)))
    (make-edge-ref e r (modulo (+ f 1) 2) id)))

(define (onext e-ref)
  (let ((e (get-edge-record e-ref))
        (r (get-rot-deg e-ref))
        (f (get-flip-deg e-ref)))
    (if (zero? f)
        ((e (modulo (+ r f) 4)) 'next)
        (flip (rot ((e (modulo (+ r f) 4)) 'next))))))


;;; Other edge functions:

(define inv-flip flip)

(define (sym e-ref)
  (rot (rot e-ref)))

(define (inv-rot e-ref)
  (rot (rot (rot e-ref))))

(define (dual e-ref)
  (rot (flip e-ref)))

(define (lnext e-ref)
  (rot (onext (inv-rot e-ref))))

(define (rnext e-ref)
  (inv-rot (onext (rot e-ref))))

(define (dnext e-ref)
  (sym (onext (sym e-ref))))

(define (oprev e-ref)
  (rot (onext (rot e-ref))))

(define (lprev e-ref)
  (sym (onext e-ref)))

(define (rprev e-ref)
  (onext (sym e-ref)))

(define (dprev e-ref)
  (inv-rot (onext (inv-rot e-ref))))


;;; Vertices are represented by an out-going edge, and faces are represented by
;;; an out-going edge in the dual diagram:

(define (org-ring e-ref)
  e-ref)

(define (left-ring e-ref)
  (org-ring (inv-rot e-ref)))

(define (right-ring e-ref)
  (org-ring (rot e-ref)))

(define (dest-ring e-ref)
  (org-ring (sym e-ref)))


;;; Basic topological operators:

(define (make-edge)
  (let* ((id (dynamic-table-size *delaunay-edges*))
         (e (make-edge-ref (make-edge-record id) 0 0 id)))
    (dynamic-table-add *delaunay-edges* e)
    e))

(define (splice a b)
  (if (not (equal? b (flip (onext a))))
      (let* ((alpha (rot (onext a)))
             (beta (rot (onext b)))

             (a-onext (onext a))
             (b-onext (onext b))
             (alpha-onext (onext alpha))
             (beta-onext (onext beta)))

        (set-onext! a b-onext beta)
        (set-onext! b a-onext alpha)
        (set-onext! alpha beta-onext b)
        (set-onext! beta alpha-onext a))))

(define (set-onext! e-ref new alt)
  (let ((f (get-flip-deg e-ref)))
    (if (zero? f)
        (let ((e (get-edge-record e-ref))
              (r (get-rot-deg e-ref)))
          (((e (modulo (+ r f) 4)) 'set-next!) new))
        (let* ((e-ref (rot (flip e-ref)))
               (e (get-edge-record e-ref))
               (r (get-rot-deg e-ref))
               (f (get-flip-deg e-ref)))
          (((e (modulo (+ r f) 4)) 'set-next!) (flip alt))))))


;;; List all edges:

(define (list-edges)
  (if (list? *delaunay-edges*)
      *delaunay-edges*
      (let ((size (dynamic-table-size *delaunay-edges*)))
        (let loop ((l '()) (i 0))
          (if (< i size)
            (let ((e (dynamic-table-fetch *delaunay-edges* i)))
              (if (eq? e 'deleted)
                (loop l (+ i 1))
                (loop (cons e l) (+ i 1))))
            (begin
             (set! *delaunay-edges* l)
             l))))))
