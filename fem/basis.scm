;;; We need to provide a common structure for handling basis functions.  This
;;; way, all that the user needs to change in order to change basis functions
;;; is the constructor passed into the domain constructor.

(declare (usual-integrations))


;;; Basis functions need to carry around their own methods:

(define (package-basis-function-methods
	 type rep eval + - * scalar*)
  (vector type rep eval + - * scalar* '()))

(define (basis-function? f)
  (and (vector? f)
       (= (vector-length f) 7)))

(define (basis:type f)
  (if (basis-function? f)
      (vector-ref f 0)
      #f))

(define (basis:get-rep f)
  (vector-ref f 1))

(define (basis-function->function f)
  (vector-ref f 2))

(define (basis:binary+ f g)
  (if (basis:same-type? f g)
      ((vector-ref f 3) g)
      (error "Cannot add basis functions of different types.")))

(define (basis:binary- f g)
  (if (basis:same-type? f g)
      ((vector-ref f 4) g)
      (error "Cannot subtract basis functions of different types.")))

(define (basis:binary* f g)
  (if (number? f)
      (if (number? g)
	  (* f g)
	  (basis:scalar* f g))
      (if (number? g)
	  (basis:scalar* g f)
	  (if (basis:same-type? f g)
	      ((vector-ref f 5) g)
	      (error "Cannot multiply basis functions of different types.")))))

(define (basis:scalar* a f)
  ((vector-ref f 6) a))

(define (basis:install-extra f tag datum)
  (let ((result (assq tag (vector-ref f 7))))
    (if result
	(set-cdr! result datum)
	(vector-set! f 7 (cons (cons tag datum) (vector-ref f 7))))))

(define (basis:get-extra f tag)
  (let ((result (assq tag (vector-ref f 7))))
    (if result
	(cdr result)
	#f)))


;;; Derived from the basic methods:

(define (basis:same-type? f g)
  (eq? (basis:type f) (basis:type g)))

(define (evaluate-basis-function f p)
  ((basis-function->function f) p))

(define (basis:+ f . rest)
  (let loop ((f f) (l rest))
    (if (null? l)
	f
	(loop (basis:binary+ f (car l)) (cdr l)))))

(define (basis:- f . rest)
  (if (null? rest)
      f
      (basis:binary- f (apply basis:+ rest))))

(define (basis:* f . rest)
  (let loop ((f f) (l rest))
    (if (null? l)
	f
	(loop (basis:binary* f (car l)) (cdr l)))))

(define (basis:dot v w)
  (let ((n (vector-length v)))
    (let loop ((i 1) (result (basis:* (vector-ref v 0) (vector-ref w 0))))
      (if (< i n)
	  (loop (+ i 1)
		(basis:+ result (basis:* (vector-ref v i) (vector-ref w i))))
	  result))))
