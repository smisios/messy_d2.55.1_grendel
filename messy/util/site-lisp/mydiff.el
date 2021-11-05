;; mydiff.el --- diff mode for GNU Emacs 20
;; (c) Rolf Sander <sander@mpch-mainz.mpg.de>
;; Time-stamp: <2005-05-25 20:10:56 sander>
 
;; to activate it copy mydiff.el to a place where emacs can find it and then
;; add "(require 'mydiff)" to your .emacs startup file
 
;; This program is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation; either version 2 of the License, or
;; (at your option) any later version.
;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.
;; You should have received a copy of the GNU General Public License
;; along with this program; if not, write to the Free Software
;; Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

;; start mydiff-mode automatically when loading a *.dif file
(setq auto-mode-alist
  (cons '("\\.dif\\'" . mydiff-mode) auto-mode-alist))

(setq mydiff-font-lock-keywords
 (list
  '("^\\+\\+\\+.*"                . font-lock-function-name-face)
  '("--> has changed everywhere!" 0 font-lock-warning-face t)
  '("--> cannot be associated unambiguously!" 0 font-lock-warning-face t)
  '("^---.*"                      . font-lock-comment-face) ; comment
  '("^<.*"                        . font-lock-keyword-face) ; 1st file
  '("^>.*"                        . font-lock-string-face)  ; 2nd file
  )
)

(defvar mydiff-mode-map () 
  "Keymap used in mydiff mode.")

(if mydiff-mode-map
    ()
  (setq mydiff-mode-map (make-sparse-keymap))
  (define-key mydiff-mode-map [f5] '(lambda () (interactive)
       (nonincremental-re-search-backward "^\\+\\+\\+ ")))
  (define-key mydiff-mode-map [f8] '(lambda () (interactive)
       (nonincremental-re-search-forward "^\\+\\+\\+ ")))
)

(defun mydiff-mode ()
  "Major mode for editing mydiff code.
Turning on mydiff mode calls the value of the variable `mydiff-mode-hook'
with no args, if that value is non-nil.

Command Table:
\\{mydiff-mode-map}"
  (interactive)
  (make-local-variable 'font-lock-defaults)
  (setq font-lock-defaults '((mydiff-font-lock-keywords) t t))
;;   (setq font-lock-defaults '((mydiff-font-lock-keywords)
;; 			     t t ((?/ . "$/"))))
  (use-local-map mydiff-mode-map)
  (setq mode-name "mydiff")
  (setq major-mode 'mydiff-mode)
  (turn-on-font-lock)
  (run-hooks 'mydiff-mode-hook)
)

(provide 'mydiff)

;; mydiff.el ends here
