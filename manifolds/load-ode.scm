;;; Load the PDE stuff separately (so that the parts of the program *not* under
;;; development can still be used without these definitions).

(let ((ode-core
       '("ode"
	 "ode-fast"
	 "lagrange"
	 "hamilton"
	 "fields"
	 "rigid-fields"
	 "ode-examples"))

      (ode-scmutils
       '("rigid")))

  (load "load-main")

  (if (not *using-scmutils?*)
      (begin
	(newline)
	(display "*** Warning: ODE code works better with ScmUtils!")
	(newline))
      (for-each load ode-scmutils))

  (for-each load ode-core))
