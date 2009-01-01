;;; Miscellaneous mathematical helpers that are useful:

(declare (usual-integrations))


;;; Some combinatorial things:

(define (choose n r)

  ;; Compute nCr:

  (/ (factorial n) (factorial r) (factorial (- n r))))

(define (slow-factorial n)
  (let loop ((n n) (result 1))
    (if (> n 1)
	(loop (- n 1) (* result n))
	result)))

(define factorial (simple-memoize slow-factorial 100))

(define (pairs l)
  (let loop ((l l) (result '()))
    (if (null? l)
	result
	(loop (cdr l)
	      (let ((a (car l)))
		(let loop ((l (cdr l)) (result result))
		  (if (null? l)
		      result
		      (loop (cdr l) (cons (list a (car l)) result)))))))))


;;; Forming the list of all sublists of L of length N is a bit more
;;; complicated.

(define (choose-sublists l n)
  (if (or (null? l) (<= n 0))
      '(())
      (let loop ((l l) (n n) (k (- (length l) n)) (result '()))
	(cond ((< k 0) result)
	      ((zero? k) (cons l result))
	      ((= n 1) (append (map list l) result))
	      (else
	       (let ((first (car l))
		     (rest (cdr l)))
		 (append
		  result
		  (loop rest n (- k 1)
			(map (lambda (sublist) (cons first sublist))
			     (loop rest (- n 1) k '()))))))))))


;;; This procedure converts references to entries in an NxN symmetric matrix
;;; into a vector representation.

(define (symmetric->vector-index i j)
  (if (>= i j)
      (if (= i 0)
	  0
	  (+ (choose (+ i 1) 2) j))
      (if (= j 0)
	  0
	  (+ (choose (+ j 1) 2) i))))


;;; Why isn't this built into Scheme?

(define (all-but list item)
  (let loop ((head '()) (tail list))
    (if (null? tail)
	list
	(if (eq? (car tail) item)
	    (append (reverse head) (cdr tail))
	    (loop (cons (car tail) head) (cdr tail))))))


;;; It's useful to find the bounding box of a finite subset of an Euclidean
;;; space:

(define (bounding-box nodes get-coords)
  (let* ((p (get-coords (car nodes)))
	 (dim (vector-length p))
	 (best (let ((l (vector->list p)))
		 (list->vector (map cons l l)))))

    (let loop ((nodes (cdr nodes)))
      (if (null? nodes)
	  (let ((l (vector->list best)))
	    (append (map car l) (map cdr l)))
	  (let ((p (get-coords (car nodes))))
	    (do ((i 0 (+ i 1)))
		((>= i dim))
	      (let ((pair (vector-ref best i))
		    (val (vector-ref p i)))
		(cond ((< val (car pair)) (set-car! pair val))
		      ((> val (cdr pair)) (set-cdr! pair val)))))
	    (loop (cdr nodes)))))))
