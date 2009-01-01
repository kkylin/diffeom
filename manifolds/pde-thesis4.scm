;;; In fixing a previous bug, I introduced another bug.

(load "pde-test")


;;; Here goes again:

(test-1 '(rectangular 10 5) '(spherical 5 10) 10000 1.9 "Data/thesis/test1a")
(test-2 '(rectangular 10 5) '(spherical 5 10) 10000 1.9 "Data/thesis/test2a")

(test-1 '(rectangular 14 7) '(spherical 7 14) 11000 1.9 "Data/thesis/test1b")
(test-2 '(rectangular 14 7) '(spherical 7 14) 11000 1.9 "Data/thesis/test2b")

(test-1 '(rectangular 18 9) '(spherical 9 18) 12000 1.9 "Data/thesis/test1c")
(test-2 '(rectangular 18 9) '(spherical 9 18) 12000 1.9 "Data/thesis/test2c")

(test-1 '(rectangular 22 11) '(spherical 11 22) 13000 1.9 "Data/thesis/test1d")
(test-2 '(rectangular 22 11) '(spherical 11 22) 13000 1.9 "Data/thesis/test2d")

(test-1 '(rectangular 26 13) '(spherical 13 26) 14000 1.9 "Data/thesis/test1e")
(test-2 '(rectangular 26 13) '(spherical 13 26) 14000 1.9 "Data/thesis/test2e")

(test-1 '(rectangular 30 15) '(spherical 15 30) 15000 1.9 "Data/thesis/test1f")
(test-2 '(rectangular 30 15) '(spherical 15 30) 15000 1.9 "Data/thesis/test2f")

(test-1 '(rectangular 34 17) '(spherical 17 34) 16000 1.9 "Data/thesis/test1g")
(test-2 '(rectangular 34 17) '(spherical 17 34) 16000 1.9 "Data/thesis/test2g")

(test-1 '(rectangular 38 19) '(spherical 19 38) 17000 1.9 "Data/thesis/test1h")
(test-2 '(rectangular 38 19) '(spherical 19 38) 17000 1.9 "Data/thesis/test2h")

(test-1 '(rectangular 42 21) '(spherical 21 42) 18000 1.9 "Data/thesis/test1i")
(test-2 '(rectangular 42 21) '(spherical 21 42) 18000 1.9 "Data/thesis/test2i")

(test-1 '(rectangular 46 23) '(spherical 23 46) 19000 1.9 "Data/thesis/test1j")
(test-2 '(rectangular 46 23) '(spherical 23 46) 19000 1.9 "Data/thesis/test2j")

(test-1 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test1k")
(test-2 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test2k")
