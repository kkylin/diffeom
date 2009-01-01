;;; Just the last two tests, which take more memory, from pde-collect.scm.

(load "pde-test")


(define test-5
  (pde:experiment-too 'pde:make-simple-domain
		      'combine-equations-with-overlap))

(define test-6
  (pde:experiment-too 'pde:make-simple-domain
		      'combine-equations-using-cmpgrd))

(test-5 '(rectangular 60 30) '(spherical 37 61) 10000 1.7 "Data/disc/test5d")
(test-6 '(rectangular 60 30) '(spherical 37 61) 10000 1.7 "Data/disc/test6d")
