;;; Just the last two tests, which take more memory, from pde-collect.scm.

(load "pde-test")


;;; Define the test procedures:

(define test-7
  (pde:experiment-too 'pde:make-simple-domain
		      'combine-equations-with-overlap1))

(define test-8
  (pde:experiment-too 'pde:make-simple-domain
		      'combine-equations-with-overlap2))


;;; Run the experiments (sorted by size, not test):

(test-7 '(rectangular 60 30) '(spherical 37 61) 10000 1.7 "Data/disc/test7d")
(test-8 '(rectangular 60 30) '(spherical 37 61) 10000 1.7 "Data/disc/test8d")
