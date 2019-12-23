G3Sim_notes2.pdf: doc/G3Sim_notes2.md
	pandoc doc/G3Sim_notes2.md -t beamer -o doc/G3Sim_notes2.pdf
doc/G3Sim_notes2.html: doc/G3Sim_notes2.md
	pandoc --self-contained doc/G3Sim_notes2.md -o doc/G3Sim_notes2.html
