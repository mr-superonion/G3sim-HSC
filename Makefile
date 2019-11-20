doc/HSM_Calib_notes.pdf: doc/HSM_Calib_notes.md 
	pandoc doc/HSM_Calib_notes.md -t beamer -o doc/HSM_Calib_notes.pdf 
doc/HSM_Calib_notes.html: doc/HSM_Calib_notes.md 
	pandoc -f markdown -s --webtex --highlight-style pygments doc/HSM_Calib_notes.md -t html -o doc/HSM_Calib_notes.html
