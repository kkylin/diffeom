;;; This file defines (again) some useful vector algebra procedures.  And
;;; SQUARE is *always* useful...

(declare (usual-integrations))


;;; Vector operations.  For completeness, we provide procedures that operate on
;;; objects of the same *shape* (made up of lists of lists of vectors, etc.),
;;; which *may* be useful for working with product structures.

(define (generalize-binary-operation combine-vectors
				     combine-numbers
				     combine-structs
				     null
				     error-string)
  (let ((report-error
	 (let ((string (string-append "Objects have different sizes or are of"
				      " the wrong type. -- "
				      error-string)))
	   (lambda ()
	     (error string))))

	(generalized-op
	 (lambda (v w)
	   (if (list? v)
	       (if (list? w)
		   (if (null? v)
		       (if (null? w)
			   null
			   (report-error))
		       (if (null? w)
			   (report-error)
			   (combine-structs (generalized-op (car v) (car w))
					    (generalized-op (cdr v) (cdr w)))))
		   (report-error))
	       (if (vector? v)
		   (if (vector? w)
		       (combine-vectors v w)
		       (report-error))
		   (if (number? v)
		       (if (number? w)
			   (combine-numbers v w)
			   (report-error))
		       (report-error)))))))
    generalized-op))


;;; Vector addition:

(define (vector:binary+ v1 v2)
  (if (= (vector-length v1) (vector-length v2))
      (let* ((len (vector-length v1))
	     (v (make-vector len)))
	(let loop ((i 0))
	  (if (< i len)
	      (begin
		(vector-set! v i (+ (vector-ref v1 i) (vector-ref v2 i)))
		(loop (+ i 1)))
	      v)))
      (error "Cannot add vectors of different dimensions. -- VECTOR:+")))

(define (vector:+ v . vlist)
  (let loop ((v v) (vlist vlist))
    (if (null? vlist)
	v
	(loop (vector:binary+ v (car vlist)) (cdr vlist)))))

(define vector:general+
  (generalize-binary-operation vector:binary+ + cons '() "VECTOR:+"))


;;; Vector subtraction:

(define (vector:binary- v1 v2)
  (if (= (vector-length v1) (vector-length v2))
      (let* ((len (vector-length v1))
	     (v (make-vector len)))
	(let loop ((i 0))
	  (if (< i len)
	      (begin
		(vector-set! v i (- (vector-ref v1 i) (vector-ref v2 i)))
		(loop (+ i 1)))
	      v)))
      (error "Cannot subtract vectors of different lengths. -- VECTOR:-")))

(define (vector:- v . vlist)
  (if (null? vlist)
      (vector:* -1 v)
      (let loop ((v v) (vlist vlist))
	(if (null? vlist)
	    v
	    (loop (vector:binary- v (car vlist)) (cdr vlist))))))

(define vector:general-
  (generalize-binary-operation vector:binary- - cons '() "VECTOR:-"))

;;; Scalar multiplication:

(define (vector:* a v)
  (let* ((len (vector-length v))
	 (w (make-vector len)))
    (let loop ((i 0))
      (if (< i len)
	  (begin
	    (vector-set! w i (* a (vector-ref v i)))
	    (loop (+ i 1)))
	  w))))

(define (vector:general* a v)
  (if (list? v)
      (if (null? v)
	  '()
	  (cons (vector:general* a (car v)) (vector:general* a (cdr v))))
      (if (vector? v)
	  (vector:* a v)
	  (if (number? v)
	      (* a v)
	      (error "Object is not a vector! -- VECTOR:*")))))


;;; Standard euclidean structures:

(define (vector:dot v w)
  (let ((len (vector-length v)))
    (if (not (= len (vector-length w)))
	(error "Vectors do not have the same dimension. -- VECTOR:DOT")
	(let loop ((i 0) (sum 0.))
	  (if (< i len)
	      (loop (+ i 1) (+ sum (* (vector-ref v i) (vector-ref w i))))
	      sum)))))

(define vector:general-dot
  (generalize-binary-operation vector:dot * + 0 "VECTOR:DOT"))


;;; Solve linear systems of equations:

(define (matrix:solve-linear-system A b)
  (let* ((m (matrix-row-count A))
	 (n (matrix-column-count A))
	 (mat (make-matrix m (+ n 1))))

    (do ((i 0 (+ i 1)))
	((>= i m))
      (do ((j 0 (+ j 1)))
	  ((>= j n))
	(matrix-set! mat i j (matrix-ref A i j))))

    (do ((i 0 (+ i 1)))
	((>= i m))
      (matrix-set! mat i n (vector-ref b i)))

    (rref mat)

    (let ((result (make-vector n)))
      (do ((j 0 (+ j 1)))
	  ((>= j n))
	(vector-set! result j (matrix-ref mat j n)))
      result)))
