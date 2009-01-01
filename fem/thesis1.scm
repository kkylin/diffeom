(load "collect")

;;; Test cases (the first run somehow stopped in the middle, probably because
;;; of errors in LU-solve):

(test-3 '(15 15 2.7) 30000 1.5 "Data/thesis/test9m")
(test-3 '(15 15 2.8) 30000 1.5 "Data/thesis/test10m")
(test-3 '(15 15 2.9) 30000 1.5 "Data/thesis/test11m")
(test-3 '(15 15 3.0) 30000 1.5 "Data/thesis/test12m")
