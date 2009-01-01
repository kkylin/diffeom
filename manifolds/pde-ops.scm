;;; Let's define some differential operators so we have something to test.

(declare (usual-integrations))


;;; Differential operators (ours only act on scalar functions, not sections of
;;; vector bundles):

(define (make-operator M make-local-form)
  (vector M make-local-form #f '()))

(define (operator:get-manifold L)
  (vector-ref L 0))

;;; Some extra structures for working with FEM.  "Context" is a kluge that lets
;;; operators learn about the particular element they are working in, etc.

(define (operator:get-local-form operator)
  (let ((context (operator:get-context operator)))
    (if context
	(apply (vector-ref operator 1) context)
	#f)))

(define (operator:set-context! operator . contextual-data)
  (vector-set! operator 2 contextual-data))

(define (operator:get-context operator)
  (vector-ref operator 2))

;;; This might come in handy:

(define (operator:install-extra L tag datum)
  (let ((result (assq tag (vector-ref L 4))))
    (if result
	(set-cdr! result datum)
	(vector-set! L 4 (cons (cons tag datum) (vector-ref L 4))))))

(define (operator:get-extra L tag)
  (let ((result (assq tag (vector-ref L 4))))
    (if result
	(cdr result)
	#f)))
