;;; Memoization is frequently useful:

(declare (usual-integrations))

;;; ALWAYS useful:

(define (square z)
  (* z z))

(define (cube z)
  (* z z z))


;;; Even more useful:

(define (simple-memoize proc size)
  ;; Memoize a function whose argument is a non-negative integer:
  (let ((cache (make-vector size 'undefined)))
    (lambda (n)
      (if (>= n size)
	  (proc n)
	  (let ((val (vector-ref cache n)))
	    (if (eq? val 'undefined)
		(let ((val (proc n)))
		  (vector-set! cache n val)
		  val)
		val))))))

(define almost-zero?
  (let ((*tolerance* 1e-14))
    (lambda (z)
      (< (magnitude z) *tolerance*))))


;;; Why is this not defined elsewhere?

(define (compose f . rest)
  (let loop ((f f) (l rest))
    (if (null? l)
	f
	(loop (binary-compose f (car l)) (cdr l)))))

(define (binary-compose f g)
  (lambda (first . rest)
    (f (apply g (cons first rest)))))
