(load "load-pde")


;;; Construct a domain:

(define make-test-domain
  (pde:make-domain-without-overlaps
   disc
   make-vertices
   make-no-extra-nodes
   planar-triangulate
   '(rectangular 20 10)
   '(spherical 10 20)))


;;; Construct the elements:

(make-test-domain
 imbedded-poly-laplacian
 make-triangular-imbedded-integrator
 pde:make-imbedded-poly-basis-function)

;;; Make a matrix!

(define mat (combine-equations-without-overlap disc 0-function test-function))


;;; Print some stuff out to file:

(if #f
    (begin
      (write-line '(writing matrix to file...))

      (let ((port (open-output-file "mat")))
	(print-matrix mat port)
	(close-output-port port))))

(write-line '(done!))

(if #f
    (let ((port (open-output-file "err")))

      (write-line `(max err = ,(max-error disc x-coord-map v)) port)
      (write-line `(min err = ,(min-error disc x-coord-map v)) port)
      (write-line `(avg err = ,(avg-error disc x-coord-map v)) port)

      (newline port)
      (write-line '(computed actual) port)

      (for-each

       (lambda (node)
	 (let ((index (node:get-id node)))
	   (if (number? index)
	       (write-line `(,(vector-ref v index)
			     ,(x-coord-map node))
			   port))))

       (sort (append-map (lambda (node)
			   (if (number? (node:get-id node))
			       (list node)
			       '()))
			 (manifold:get-nodes disc))
	     (lambda (n1 n2)
	       (< (node:get-id n1) (node:get-id n2)))))

      (close-output-port port)))
