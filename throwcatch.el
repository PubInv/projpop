(require 'cl-lib)

;; Here I create the basic parameters.
;; In KMS system...
;; 1 kilogram
(setq A 0.3)

;; kilograms
(setq B 0.1)
 ;; 1 kilogram
(setq C 1)


;; Apparently S is really irrelevant to all our math!
;; meters / s
(setq S 100)



;; Here I create the recurrence relations as matehmatical funcions
;; Note: There is an omnipresent factor of "S" in each of these results
;; which I will not represent until the final step.
;; The range of these functions start at 1 -- that is probably
;; not too elegant, but that is how I set that up.

;; Note: Surprisingly, the number of bounces is completely independent of
;; the initial velocity!!! In a sense, it depends only on relative masses of the systems!

;; Would be more elegant if A, B, C where parameters (and the table!)
;; Have we got any hope of a closed form solution of the number of bounces?


;; So this is rather interesting...
;; We might now ask the quesiton: How do we obtain a closed form solution?
;; Or perhaps we should just build a few graphs...
;; Note: if we were serious about using this for propulsion, we would think
;; about the time it takes for each trip, and acceleration over time, etc.

;; Note: Increase C to B ratio by 10 tends to triple the number of trips.
;; Note: Decrease the mass ratio significantly effects final velocity!

;; I think we could produce analytics bounds and also compute the time
;; to take the trips and the distance travelled at the final collision,
;; this would be an interesting paper.
;; If we could relate this to the "rocket equation" it could be interesting.
;; What if a rocket intentionally outputs a small cup and then keeps this bouncing?
;; Is that more efficient than driving back the cup with some bigger engine?
;; The beauty is that your are reusing the propellant.
;; Probably not advantageous unless this lets you build a better motor in some way.
;; Probably more interesting is the anchored case.

;; But what is worth "writing up" : Fact that it is independent of velocity is interesting.
;; The graph itself may be interesting.
;; Really have to understand how it compares to shooting a ball at a certain speed.
;; That is actually interesting. If a rocket begins in space with a certain mass,
;; and it can shoot a ball, can it be more efficient shooting a ball this way
;; than just shooting the balls in the first place?  Solve the math for that.

(defun va (n tab)
  (if (= n 1)
      0
    (if (cl-evenp n)
	(vat (- n 1) tab)
      (/ (+
	  (* 2 B (vbt (- n 1) tab))
	  (* (- A B) (vat (- n 1) tab)))
	 (+ A B)))
    ))

(defun vb (n tab)
  (if (= n 1)
      1  ;; This would be "S"
    (if (cl-evenp n)
      (/ 
       (+ 
        (* (- B C) (vbt (- n 1) tab))
	(* 2 C (vct (- n 1) tab))
	)
       (+ B C))
      (/ 
       (+ 
        (* (- B A) (vbt (- n 1) tab))
	(* 2 A (vat (- n 1) tab))
	)
       (+ A B))
      )))


(defun vc (n tab)
  (if (= n 1)
      0
    (if (cl-oddp n)
	(vct (- n 1) tab)
      (/ (+
	  (* 2 B (vbt (- n 1) tab))
	  (* (- C B) (vct (- n 1) tab))
	    )
	 (+ B C)))
    ))

(defun will-collidep (va vb vc)
  "true if and only if B will strike A or C (assuming between them)"
  (or (< vb va) (> vb vc))
 )


(defun kinetic-energy (va vb vc)
  "return the total kinetic energy of these velocitys"
  (/
   (+
    (* A (* va va))
    (* B (* vb vb))
    (* C (* vc vc))
   2)
  ))

;; Compute a table.

(defun vat (n tab)
  "Look up the value of va from a table"
  (let ((row (assoc n tab)))
    (nth 1 row)))

(defun vbt (n tab)
  "Look up the value of vb from a table"
  (let ((row (assoc n tab)))
    (nth 2 row)))

(defun vct (n tab)
  "Look up the value of vc from a table"
  (let ((row (assoc n tab)))
    (nth 3 row)))


(defun travel-time (n apos cpos va vb vc)
  "travel-time for the nth collision"
  (if (cl-oddp n)
      ;; time to travel fro a to c, given that
      ;; c is traveling at vc
      (/ (+ (- apos) cpos)
	 (* S (- vb vc)))
    (/ (+ (- apos) cpos)
       (* S (- va vb)))
      ))


(defun table (limit)
  "compute up to limit collisions"
  (let* ((value nil)
	 (i 0)
	 (cp t)
	 (cpos 1.0)
	 (apos -1.0)
	 (napos nil)
	 (ncpos nil)
	 (travel nil)
	 (history 0.0))
    (while (and cp (< i limit))
      (let* ((n (+ i 1))
	     (ax (va n value))
	     (bx (vb n value))
	     (cx (vc n value))
	     )
	(setf tt (travel-time n apos cpos ax bx cx))
	(setf napos (+ apos (* tt ax)))
	(setf ncpos (+ cpos (* tt cx)))
	(setf cp (will-collidep ax bx cx))
	    (setq value
		  (cons
		   (list n ax bx cx
			 cp
			 (kinetic-energy ax bx cx)
			 (list
			  history
			  tt
			  (list
			 apos
			 cpos)))
		   
		   value))
      (setf i (+ i 1))
      (setf cpos ncpos)
      (setf apos napos)
      (setf history (+ history tt))
      ))
    value))

(defun format-row-list (row)
  (list
   (car row)
   (* S (nth 1 row))
   (* S (nth 2 row))
   (* S (nth 3 row))
   (nth 6 row)
   ))

(defun format-row (row)
  (let ((poss (nth 2 (nth 6 row))))
    (concat
     (format "col:   %d " (car row))
     (format "vA %.3f m/s " (* S (nth 1 row)))
     (format "vB %.3f m/s " (* S (nth 2 row)))
     (format "vC %.3f m/s " (* S (nth 3 row)))
     (format "\n")
     (format "A: %.3f m " (nth 0 poss))
     (format "C: %.3f m " (nth 1 poss))
     (format "\n")     
    )
     ))

(defun pptable (limit)
  "pretty print a table"
  (let ((tab (table limit))
	(value nil))
    (dolist (x tab value)
      (setf value (concat (format-row x) value))
      )))


(defun tsiolkovsky (ve m0 mf)
  (* ve (log (/ m0 mf))))

(defun rocket-speed ()
  (tsiolkovsky S (+ A B C) C))
