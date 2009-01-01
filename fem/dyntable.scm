(declare (usual-integrations))

(define (make-dynamic-table . argl)
  (if (or (null? argl)
	  (> (length argl) 0))
      (vector 0 (make-vector 8))
      (vector 0 (make-vector (car argl)))))

(define (dynamic-table-size table)
  (vector-ref table 0))

(define (dynamic-table-add table new-element)
  (let ((size (dynamic-table-size table))
	(real-size (vector-length (vector-ref table 1))))
    (if (= size real-size)
	(let ((old-table (vector-ref table 1))
	      (new-table (make-vector (* 2 real-size))))
	  (do ((i 0 (+ i 1)))
	      ((>= i real-size))
	    (vector-set! new-table i (vector-ref old-table i)))
	  (vector-set! table 1 new-table)))

    (vector-set! (vector-ref table 1) size new-element)
    (vector-set! table 0 (+ size 1))))

(define (dynamic-table-fetch table i)
  (if (>= i (dynamic-table-size table))
      (error "Access out of bound -- DYNAMIC-TABLE-FETCH")
      (vector-ref (vector-ref table 1) i)))

(define (dynamic-table->list table)
  (let loop ((l '()) (n (- (dynamic-table-size table) 1)))
    (if (< n 0)
	l
	(loop (cons (dynamic-table-fetch table n) l) (- n 1)))))

(define (dtableq table element)
  (let loop ((n (- (dynamic-table-size table) 1)))
    (if (< n 0)
	#f
	(if (eq? (dynamic-table-fetch table n) element)
	    n
	    (loop (- n 1))))))

(define (dynamic-table-set! table i val)
  (if (>= i (dynamic-table-size table))
      (error "Access out of bound -- DYNAMIC-TABLE-SET!")
      (vector-set! (vector-ref table 1) i val)))
