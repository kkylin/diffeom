;;; See pde-test.scm.old for more information about the errors associated with
;;; different variants of our method.


;;; Let's see how the error scales with the number of nodes:

(define (pde:experiment domain-maker combine-equations)

  ;; DOMAIN-MAKER and COMBINE-EQUATIONS should be *symbols*, not procedures.
 
  (lambda (rectangular spherical sor-steps sor-coeff port)

    ;; First, reload everything (to clear residual states in data structures).

    (load "load-pde")
    (if port (newline port))

    ;; Let's start:

    (let ((make-test-domain
	   ((evaluate-symbol domain-maker) disc
					   make-vertices
					   make-no-extra-nodes
					   planar-triangulate
					   rectangular
					   spherical))
	  (f test-function))

      (make-test-domain
       imbedded-poly-laplacian
       make-triangular-imbedded-integrator
       pde:make-imbedded-poly-basis-function)

      (let ((mat ((evaluate-symbol combine-equations)
		  disc 0-function test-function)))

	(write-line '(creating normal equations...))
	(set! mat (show-time (lambda () (sparse-normal-equations mat))))

	(write-line '(solving normal equations...))

	(let ((v (show-time (lambda () (sor mat sor-steps sor-coeff))))
	      (write (lambda (stuff)
		       (write-line stuff)
		       (if port
			   (write-line stuff port)))))
	  (write `(domain-maker = ,domain-maker))
	  (write `(combine-equations = ,combine-equations))
	  (write `(,(length (manifold:get-nodes disc)) nodes))
	  (write `(max absolute error = ,(max-error disc f v)))
	  (write `(min absolute error = ,(min-error disc f v)))
	  (write `(average absolute error = ,(avg-error disc f v)))
	  (write `(max relative error = ,(max-relative-error disc f v)))
	  (write `(min relative error = ,(min-relative-error disc f v))))))))


;;; Use #f for FILE-NAME if standard output is the only desired output port.
;;; Otherwise, the output is sent to both standard output and the named file.

(define (run-test-case test-case file-name)
  (let ((port (if file-name
		  (open-output-file file-name)
		  #f)))

    (write-line `(test case: ,test-case))
    (if port (write-line `(test case: ,test-case) port))

    (case test-case
      ((1)
       (let ((try-pde (pde:experiment 'pde:make-domain-with-overlaps
				      'combine-equations-with-overlap)))
	 (try-pde '(rectangular 10 5) '(spherical 5 10) 10000 1.9 port)
	 (try-pde '(rectangular 20 10) '(spherical 10 20) 20000 1.9 port)
	 (try-pde '(rectangular 40 25) '(spherical 25 40) 30000 1.9 port)))
      ((2)
       (let ((try-pde (pde:experiment 'pde:make-domain-with-overlaps
				      'combine-equations-with-overlap)))
	 (try-pde '(rectangular 20 10) '(spherical 10 20) 10000 1.9 port)))
      ((3)
       (let ((try-pde (pde:experiment 'pde:make-simple-domain
				      'combine-equations-using-CMPGRD)))
	 (try-pde '(rectangular 20 10) '(spherical 10 20) 10000 1.9 port)))

      ((4)
       (let ((try-pde (pde:experiment 'pde:make-simple-domain
				      'combine-equations-using-CMPGRD)))
	 (try-pde '(rectangular 10 5) '(spherical 5 10) 10000 1.9 port)
	 (try-pde '(rectangular 20 10) '(spherical 10 20) 20000 1.9 port)
	 (try-pde '(rectangular 40 25) '(spherical 25 40) 30000 1.9 port)))

      ((5)
       (let ((try-pde (pde:experiment 'pde:make-domain-with-overlaps
				      'combine-equations-with-overlap)))
	 (try-pde '(rectangular 10 5) '(spherical 5 10) 0 0 #f)))

      ((6)
       (let ((try-pde (pde:experiment 'pde:make-domain-with-small-overlaps
				      'combine-equations-without-overlaps)))
	 (try-pde '(rectangular 20 10) '(spherical 10 20) 10000 1.9 #f)))

      (else #f))

    (if port (close-output-port port))))


;;; A slight variant used for collecting data for thesis work:

(define (pde:experiment-too domain-maker combine-equations)

  ;; DOMAIN-MAKER and COMBINE-EQUATIONS should be *symbols*, not procedures.

  (lambda (rectangular spherical sor-steps sor-coeff file)

    ;; First, reload everything (to clear residual states in data structures).

    (load "load-pde")

    ;; Let's start:

    (show-time
     (lambda ()
       (let ((make-domain
	      ((evaluate-symbol domain-maker) disc
					      make-vertices
					      make-no-extra-nodes
					      planar-triangulate
					      rectangular
					      spherical))
	     (f test-function))

	 (write-line '(constructing elements...))

	 (make-domain
	  imbedded-poly-laplacian
	  make-triangular-imbedded-integrator
	  pde:make-imbedded-poly-basis-function)

	 (let ((mat ((evaluate-symbol combine-equations) disc 0-function f)))

	   (if (not (= (+ (sparse-matrix-row-count mat) 1)
		       (sparse-matrix-column-count mat)))
	       (begin
		 (write-line '(need normal equations.))
		 (set! mat (show-time
			    (lambda ()
			      (sparse-normal-equations mat))))))

	   (write-line '(solving equations...))

	   (let ((v (show-time (lambda () (sor mat sor-steps sor-coeff)))))
	     (write-line '(preparing to save data...))
	     (let ((states (node-states disc f v))
		   (port (open-output-file file)))
	       (write-line '(saving...))
	       (print-matrix states port)
	       (close-output-port port)))))))))


;;; Tests for thesis-related data:

(define test-1
  (pde:experiment-too 'pde:make-domain-without-overlaps
		      'combine-equations-without-overlap))

(define test-2
  (pde:experiment-too 'pde:make-domain-with-small-overlaps
		      'combine-equations-without-overlap))

(define test-3
  (pde:experiment-too 'pde:make-domain-with-overlaps
		      'combine-equations-with-overlap))

(define test-4
  (pde:experiment-too 'pde:make-domain-with-larger-overlaps
		      'combine-equations-with-overlap))

(define test-5
  (pde:experiment-too 'pde:make-simple-domain
		      'combine-equations-with-overlap))

(define test-6
  (pde:experiment-too 'pde:make-simple-domain
		      'combine-equations-using-cmpgrd))

(define test-7
  (pde:experiment-too 'pde:make-simple-domain
		      'combine-equations-with-overlap1))

(define test-8
  (pde:experiment-too 'pde:make-simple-domain
		      'combine-equations-with-overlap2))
