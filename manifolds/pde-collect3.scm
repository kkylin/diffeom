;;; Even more tests:

(load "pde-test")


;;; Define the test procedures:

(define test-7
  (pde:experiment-too 'pde:make-simple-domain
		      'combine-equations-with-overlap1))

(define test-8
  (pde:experiment-too 'pde:make-simple-domain
		      'combine-equations-with-overlap2))


;;; Run the experiments (sorted by size, not test):

(test-7 '(rectangular 40 15) '(spherical 20 30) 10000 1.7 "Data/disc/test7c")
(test-8 '(rectangular 40 15) '(spherical 20 30) 10000 1.7 "Data/disc/test8c")
