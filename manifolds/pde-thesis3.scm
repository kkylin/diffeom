;;; Turns out one of the bug fixes introduced an error.  The new code with the
;;; error got reloaded into the system, and...

(load "pde-test")


(test-5 '(rectangular 46 23) '(spherical 23 46) 19000 1.9 "Data/thesis/test5j")

(test-1 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test1k")
(test-2 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test2k")
(test-3 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test3k")
(test-4 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test4k")
(test-5 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test5k")
