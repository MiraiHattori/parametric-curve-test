(progn
    (setq G #2f((4 1) (1 2))) ;; eval-weight-matrix
    (setq initial-state #f(0 0)) ;; initial-state
    (print "Be careful. ce0 signature is reversed when input it into eiquadprog.")
    (setq ce0 #f(1)) ;; equality-vector
    (setq CE #2f((1)
                 (1))) ;; equality-matrix
    (setq CI #2f((1)
                 (1))) ;; inequality-matrix
    (setq ci0-min #f(-inf)) ;; inequality-min-vector
    (setq ci0-max #f(inf)) ;; inequality-max-vector
    (load "package://eus_qp/euslisp/eiquadprog.l")
    (print (solve-eiquadprog :initial-state initial-state
                             :eval-weight-matrix G
                             :equality-vector ce0
                             :equality-matrix CE
                             :inequality-matrix CI
                             :inequality-min-vector ci0-min
                             :inequality-max-vector ci0-max
                             :eqthre 1e-5)))
