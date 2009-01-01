;;; Use this file to collect data for theis work.  Based on pde-test.scm.

(load "pde-test")


;;; Define the test procedures:

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


;;; Run the experiments (sorted by size, not test):

(test-1 '(rectangular 10 5) '(spherical 5 10) 10000 1.7 "Data/disc/test1a")
(test-2 '(rectangular 10 5) '(spherical 5 10) 10000 1.7 "Data/disc/test2a")
(test-3 '(rectangular 10 5) '(spherical 5 10) 10000 1.7 "Data/disc/test3a")
(test-4 '(rectangular 10 5) '(spherical 5 10) 10000 1.7 "Data/disc/test4a")
(test-5 '(rectangular 10 5) '(spherical 5 10) 10000 1.7 "Data/disc/test5a")
(test-6 '(rectangular 10 5) '(spherical 5 10) 10000 1.7 "Data/disc/test6a")

(test-1 '(rectangular 20 10) '(spherical 10 20) 10000 1.7 "Data/disc/test1b")
(test-2 '(rectangular 20 10) '(spherical 10 20) 10000 1.7 "Data/disc/test2b")
(test-3 '(rectangular 20 10) '(spherical 10 20) 10000 1.7 "Data/disc/test3b")
(test-4 '(rectangular 20 10) '(spherical 10 20) 10000 1.7 "Data/disc/test4b")
(test-5 '(rectangular 20 10) '(spherical 10 20) 10000 1.7 "Data/disc/test5b")
(test-6 '(rectangular 20 10) '(spherical 10 20) 10000 1.7 "Data/disc/test6b")

(test-1 '(rectangular 40 15) '(spherical 20 30) 10000 1.7 "Data/disc/test1c")
(test-2 '(rectangular 40 15) '(spherical 20 30) 10000 1.7 "Data/disc/test2c")
(test-3 '(rectangular 40 15) '(spherical 20 30) 10000 1.7 "Data/disc/test3c")
(test-4 '(rectangular 40 15) '(spherical 20 30) 10000 1.7 "Data/disc/test4c")
(test-5 '(rectangular 40 15) '(spherical 20 30) 10000 1.7 "Data/disc/test5c")
(test-6 '(rectangular 40 15) '(spherical 20 30) 10000 1.7 "Data/disc/test6c")

(test-1 '(rectangular 60 30) '(spherical 37 61) 10000 1.7 "Data/disc/test1d")
(test-2 '(rectangular 60 30) '(spherical 37 61) 10000 1.7 "Data/disc/test2d")
(test-3 '(rectangular 60 30) '(spherical 37 61) 10000 1.7 "Data/disc/test3d")
(test-4 '(rectangular 60 30) '(spherical 37 61) 10000 1.7 "Data/disc/test4d")
(test-5 '(rectangular 60 30) '(spherical 37 61) 10000 1.7 "Data/disc/test5d")
(test-6 '(rectangular 60 30) '(spherical 37 61) 10000 1.7 "Data/disc/test6d")
