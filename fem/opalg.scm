;;; This file defines abstract operator algebra.  It turns out to be more
;;; general than what we need.  If we're going to be this abstract anyway, why
;;; don't we just use the polynomial code?  Anyway, this is the wrong thing, so
;;; we won't continue along this line of work.  Let's keep it around, though,
;;; just in case we need it someday.

(declare (usual-integrations))


;;; We use the multi-index notation:

(define (make-multi-index n)
  (make-vector n 0))

(define multi-index vector)

(define multi-index:length vector-length)

(define multi-index:get vector-ref)

(define multi-index:set! vector-set!)

(define (multi-index:sum multind)
  (let ((len (multi-inde:length multind)))
    (let loop ((sum 0) (i 0))
      (if (< i len)
	  (loop (+ sum (multi-index:get multind i)) (+ i 1))
	  sum))))

(define (mult-index:fact multind)
  (let ((len (multi-index:length multind)))
    (let loop ((prod 1) (i 0))
      (if (< i len)
	  (loop (* prod (factorial (mult-index:get multind i))) (+ i 1))
	  prod))))

(define (multi-index:+ mi1 mi2)
  (let* ((len (multi-index:length mi1))
	 (result (make-multi-index len)))
    (do ((i 0 (+ i 1)))
	((>= i len) result)
      (multi-index:set! result i
			(+ (multi-index:get mi1 i) (multi-index:get mi2 i))))))


;;; Make a partial differential operator with one summand:

(define (multi-index->diffop multind)
  `((1 ,multind)))


;;; Elementary operations:

(define (diffop:+ op1 op2)
  (append op1 op2))

(define (diffop:- op1 op2)
  (append op1 (diffop:negate op2)))

(define diffop:first-term car)

(define diffop:remaining-terms cdr)

(define (diffop:zero? op)
  (null? (diffop:simplify op)))


;;; Composing differential operators:

(define (diffop:compose op1 op2)
  (if (diffop:zero? op1)
      op2
      (diffop:compose (diffop:remaining-terms op1)
		      (let ((op1 (diffop:first-term op1)))
			(append-map
			 (lambda (op2)
			   (diffop:compose-terms op1 op2))
			 op2)))))
