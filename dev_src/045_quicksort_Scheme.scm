;;; InPlaceQuicksort.scm
;;; Scheme — In-place quicksort on a vector.
;;;
;;; Features
;;; • Reads integers (one per line or whitespace-separated) from a file
;;;   if a path is provided as the first CLI arg; otherwise from STDIN.
;;; • In-place quicksort using an iterative stack + Hoare partition
;;;   with median-of-three pivot selection.
;;; • Prints sorted numbers to STDOUT, one per line.
;;;
;;; Tested with GNU Guile; should work on many Scheme impls that provide
;;; `command-line`. If your Scheme lacks `command-line`, call (main) or
;;; (main "path/to/file") from the REPL.

;; ----------------------------
;; Utility: CLI arg (best-effort)
;; ----------------------------
(define (first-cli-arg)
  (cond
    ;; GNU Guile
    ((defined? 'command-line)
     (let ((argv (command-line)))
       (if (> (length argv) 1) (list-ref argv 1) #f)))
    ;; Racket compatibility (when run as #!r6rs / racket -I r5rs) — not guaranteed
    ((defined? 'current-command-line-arguments)
     (let ((vec (current-command-line-arguments)))
       (if (> (vector-length vec) 0) (vector-ref vec 0) #f)))
    (else #f)))

;; ----------------------------
;; Input: read numbers into vector
;; ----------------------------
(define (read-all-ints port)
  ;; Read s-expressions; keep only exact integers or inexact that are whole.
  (let loop ((xs '()))
    (let ((x (read port)))
      (if (eof-object? x)
          (list->vector (reverse xs))
          (cond
            ((integer? x)      (loop (cons x xs)))
            ((and (number? x) (= x (inexact->exact (round x)))))
             (loop (cons (inexact->exact (round x)) xs)))
            (else (error "Non-integer token in input:" x)))))))

;; ----------------------------
;; Quicksort helpers (in place)
;; ----------------------------
(define (swap! v i j)
  (let ((tmp (vector-ref v i)))
    (vector-set! v i (vector-ref v j))
    (vector-set! v j tmp)))

(define (median3 a b c)
  ;; Return the median of three numbers.
  (cond
    ((or (and (<= a b) (<= b c)) (and (<= c b) (<= b a))) b)
    ((or (and (<= b a) (<= a c)) (and (<= c a) (<= a b))) a)
    (else c)))

(define (hoare-partition! v lo hi)
  ;; Hoare partition with median-of-three pivot.
  ;; Returns index j such that [lo..j] <= pivot and [j+1..hi] >= pivot.
  (let* ((mid (+ lo (quotient (- hi lo) 2)))
         (pivot (median3 (vector-ref v lo)
                         (vector-ref v mid)
                         (vector-ref v hi)))
         (i (- lo 1))
         (j (+ hi 1)))
    (let loop ()
      (begin
        (let nxt-i ()
          (set! i (+ i 1))
          (if (< (vector-ref v i) pivot) (nxt-i)))
        (let nxt-j ()
          (set! j (- j 1))
          (if (> (vector-ref v j) pivot) (nxt-j)))
        (if (>= i j)
            j
            (begin (swap! v i j) (loop)))))))

(define (quicksort-in-place! v)
  ;; Iterative quicksort using an explicit stack of (lo . hi) pairs.
  (let loop ((stack (list (cons 0 (- (vector-length v) 1)))))
    (if (null? stack)
        v
        (let* ((seg (car stack))
               (lo (car seg))
               (hi (cdr seg))
               (rest (cdr stack)))
          (if (< lo hi)
              (let ((p (hoare-partition! v lo hi)))
                ;; Push subranges
                (loop (let ((s rest))
                        (if (< lo p) (set! s (cons (cons lo p) s)))
                        (if (< (+ p 1) hi) (set! s (cons (cons (+ p 1) hi) s)))
                        s)))
              (loop rest))))))

;; ----------------------------
;; Output
;; ----------------------------
(define (print-ints v)
  (let ((n (vector-length v)))
    (let iter ((i 0))
      (when (< i n)
        (display (vector-ref v i)) (newline)
        (iter (+ i 1))))))

;; ----------------------------
;; Entrypoint
;; ----------------------------
(define (main . maybe-arg)
  (let* ((path (cond
                 ((and (pair? maybe-arg) (string? (car maybe-arg))) (car maybe-arg))
                 ((first-cli-arg) => (lambda (s) s))
                 (else #f)))
         (vec  (if path
                   (call-with-input-file path read-all-ints)
                   (read-all-ints (current-input-port)))))
    (quicksort-in-place! vec)
    (print-ints vec)))

;; Auto-run if loaded as a script in many Schemes
;; (harmless if your implementation ignores `command-line`).
(when (and (defined? 'command-line)
           ;; crude heuristic: if file was invoked directly, try running:
           #t)
  (when (null? (interaction-environment))
    (main)))

