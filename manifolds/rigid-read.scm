;;; Generate some plots from data!  Of course, the two runs are using different
;;; formats, which makes it kind of hard...

(declare (usual-integrations))


(define (eof? port)
  (eof-object? (peek-char port)))

(define (read-number port)
  (let loop ((c (read-char port)) (string '()))
    (if (or (eof-object? c) (memq c '(#\space #\tab #\newline)))
	(string->number (list->string (reverse string)))
	(if (memq c (string->list "-.0123456789e"))
	    (loop (read-char port) (cons c string))
	    (loop (read-char port) string)))))

(define (read-vector n port)
  (let ((v (make-vector n)))
    (do ((i 0 (+ i 1)))
	((>= i n) v)
      (vector-set! v i (read-number port)))))

(define (read-regular-state port)
  (let ((t (read-number port))
	(x (read-vector 3 port))
	(v (read-vector 3 port)))
    (->state t x v)))

(define read-manifold-state
  (let ((charts (manifold:get-finite-atlas TSO3))
	(Tchart (make-tangent-chart euler-angles)))
    (lambda (port)
      (let* ((t (read-number port))
	     (chart (read-number port))
	     (coords (read-vector 6 port))

	     (v (chart:point->coords
		 (chart:coords->point coords (list-ref charts chart))
		 Tchart))

	     (psi (vector-ref v 0))
	     (theta (vector-ref v 1))
	     (phi (vector-ref v 2))
	     (psidot (vector-ref v 3))
	     (thetadot (vector-ref v 4))
	     (phidot (vector-ref v 5)))
	(->state t (vector theta phi psi)
		 (vector thetadot phidot psidot))))))

(define (read-regular-file filename)
  (let ((port (open-input-file filename)))
    (let loop ((states '()) (count 0) (total 0))
      (if (eof? port)
	  (begin
	    (close-input-port port)
	    (sort states (lambda (x y) (< (state->t x) (state->t y)))))
	  (if (> count 100)
	      (begin
		(write-line `(read ,total states))
		(loop (cons (read-regular-state port) states) 0 (+ total 1)))
	      (loop (cons (read-regular-state port) states)
		    (+ count 1) (+ total 1)))))))

(define (read-manifold-file filename)
  (let ((port (open-input-file filename)))
    (let loop ((states '()) (count 0) (total 0))
      (if (eof? port)
	  (begin
	    (close-input-port port)
	    (sort states (lambda (x y) (< (state->t x) (state->t y)))))
	  (if (> count 100)
	      (begin
		(write-line `(read ,total states))
		(loop (cons (read-manifold-state port) states) 0 (+ total 1)))
	      (loop (cons (read-manifold-state port) states)
		    (+ count 1) (+ total 1)))))))

(define (read-pendulum-file filename)
  (let ((port (open-input-file filename)))
    (let loop ((result '()) (count 0) (total 0))
      (if (eof? port)
	  (begin
	    (close-input-port port)
	    (reverse result))
	  (if (> count 100)
	      (begin
		(write-line `(read ,total states))
		(loop (cons (read-pendulum-state port) result) 0 (+ total 1)))
	      (loop (cons (read-pendulum-state port) result)
		    (+ count 1) (+ total 1)))))))

(define (read-pendulum-state port)
  (read-number port)
  (let* ((x (read-number port))
	 (y (read-number port))
	 (z (read-number port))
	 (px (read-number port))
	 (py (read-number port))
	 (pz (read-number port)))
    (imbedding->cotangent S^2 (vector x y z) (vector px py pz))))
