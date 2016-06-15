(require 'cl-lib)

;; Here I create the basic parameters.
;; In KMS system...
;; 1 kilogram
(setq A 0.1)

;; kilograms
(setq B 0.2)
 ;; 1 kilogram
(setq C 1)


;; Apparently S is really irrelevant to all our math!
;; meters / s
(setq S 10)


;; Interesting note: There are values in which the speed
;; approaches the rocket-speed, even though we are applying
;; energy only to B, not A + B.  For example:
;; (A,B,C) = (0.1,0.2,1.0)
;; Seems to have a Vc of 33 m/s, while the rocket speed is 26.23!
;; We need to compare this to a regime of firing separately.
;; I am not sure what this means.  I certainly
;; believe that if we built a "multi-shot" sequence we
;; can outdo the rocket equation!  Perhaps that is known.
;; Perhaps it is completely revolutionary!

;; TODO: Understand if this is really true and how significant it is.
;; Question: if we imagine the first object reaceding, does this
;; still work?  Or it is too difficult to bounce?

;; I'm afraid there is now way around the fact that we need
;; a more general physical and computational model.
;; The basis of this one-dimenstional model could be to
;; having a lot of objects on the track interacting.
;; we need an function that predicts the "next" collision...
;; But computationally this is easy because you just look
;; at each particle up and down the track, so O(n).
;;  

;; PRORITY: Turn it better into a rocket problem.
;; Imagine firing A backwards, and the B, and in general
;; being able to set up a regime of masses and speeds.
;; I predict that we will want to use smaller masses at
;; higher speeds in the later operations.

;; This would mean expanding this to have a physical state
;; that represents many masses.


;; Here I create the recurrence relations as matehmatical funcions
;; Note: There is an omnipresent factor of "S" in each of these results
;; which I will not represent until the final step.
;; The range of these functions start at 1 -- that is probably
;; not too elegant, but that is how I set that up.

;; Note: Surprisingly, the number of bounces is completely independent of
;; the initial velocity!!! In a sense, it depends only on relative masses of the systems!

;; Would be more elegant if A, B, C where parameters (and the table!)
;; Have we got any hope of a closed form solution of the number of bounces?
;; Note: There is a sense in which it is optimal to cut the mass down
;; so that the you just barely recover the reaction mass.
;; It seems that an optimal situation may be a pulsed output of a
;; heavy mass that you bounce things off until you abandon a mass,
;; and then send out another mass.


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

;; Ideas for physical representation of state:
;; Keep a simple list of objects ordered by physical position.
;; Be able to compute the time of the next interaction.
;; One "turn" advances to the next interaction, updating positions and velocities
;; accordingly. An Interaction function can later be introduced to support
;; non-elastic collisions.



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
      )
    (/ (* S (nth 3 (car tab))) (rocket-speed))
    ))

(defun find-rocket-speed-max (limit)
  ;; iterate over A and B by interval seeking highest ratio over rocket speed
  (let ((i 0.0)
	(j 0.0)
	(max ()))
    (while (< j 1.0)
      (while (< i 1.0)
	(let ((ratio (pptable limit)))
	  (pp ratio)
	  (pp " ")
	  (if (or (eq max nil) (< (car max) ratio))
	      (setf max (list ratio 'A i 'B j)))
	  )
	(setf i (+ i 0.1))
	)
      (setf j (+ j 0.1))
      )
    max))



(defun tsiolkovsky (ve m0 mf)
  (* ve (log (/ m0 mf))))

(defun rocket-speed ()
  (tsiolkovsky S (+ A B C) C))
