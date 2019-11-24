doc/HSM_Calib_notes.pdf: doc/HSM_Calib_notes.md 
	pandoc doc/HSM_Calib_notes.md -t beamer -o doc/HSM_Calib_notes.pdf 
doc/HSM_Calib_notes.html: doc/HSM_Calib_notes.md 
	pandoc --self-contained doc/HSM_Calib_notes.md -o doc/HSM_Calib_notes.html
