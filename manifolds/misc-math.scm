(declare (usual-integrations))


;;; Some definitions that are always useful:

(define (square z)
  (* z z))

(define (compose f . rest)
  (let loop ((f f) (l rest))
    (if (null? l)
	f
	(loop (binary-compose f (car l)) (cdr l)))))

(define (binary-compose f g)
  (lambda (first . rest)
    (f (apply g (cons first rest)))))

(define pi (* 4 (atan 1)))
(define -pi (- pi))

(define (identity x)
  x)

(define (relative-error val ref)
  (if (zero? ref)
      val
      (/ (- val ref) ref)))
