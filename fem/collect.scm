;;; Collect some data for thesis work.  First, load the FEM programs and set
;;; the speed of light to 1.


;;; Here's how we run experiments:

(define (make-experiment domain-maker)
  (lambda (argl filename)

    ;; Reload to clear hidden states:

    (load "load")
    (set! *wave-constant* 1.)
    (let ((make-domain (evaluate-symbol domain-maker)))

      (write-line '(constructing domain...))
      (let ((nodes (show-time (lambda () (apply make-domain argl)))))
	(write-line `(,(vector-length nodes) nodes constructed))
	(write-line '(constructing matrix...))
	(let ((mat (show-time
		    (lambda ()
		      (sparse->matrix (fem 0-function nodes wave))))))
	  (write-line `(matrix size = ,(matrix-size mat)))
	  (write-line '(solving equations...))
	  (let ((v (show-time (lambda () (lu-solve mat)))))
	    (write-line '(computing results...))
	    (let ((results (compute-results nodes v wave)))
	      (write-line '(saving...))
	      (let ((port (open-output-file filename)))
		(print-matrix results port)
		(close-output-port port)))))))))

(define (make-relaxing-experiment domain-maker)
  (lambda (argl sor-count sor-factor filename)

    ;; Reload to clear hidden states:

    (load "load")
    (set! *wave-constant* 1.)
    (let ((make-domain (evaluate-symbol domain-maker)))

      (write-line '(constructing domain...))
      (let ((nodes (show-time (lambda () (apply make-domain argl)))))
	(write-line `(,(vector-length nodes) nodes constructed))
	(write-line '(constructing matrix...))
	(let ((mat (show-time
		    (lambda ()
		      (sparse-normal-equations
		       (fem 0-function nodes wave))))))
	  (write-line `(matrix size = ,(sparse-matrix-size mat)))
	  (write-line '(relaxing...))
	  (let ((v (show-time (lambda () (sor mat sor-count sor-factor)))))
	    (write-line '(computing results...))
	    (let ((results (compute-results nodes v wave)))
	      (write-line '(saving...))
	      (let ((port (open-output-file filename)))
		(print-matrix results port)
		(close-output-port port)))))))))


;;; A couple of tests:

(define test-1 (make-experiment 'make-true-hat-domain))
(define test-2 (make-experiment 'make-bent-domain))
(define test-3 (make-relaxing-experiment 'make-true-hat-domain))
