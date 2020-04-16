OCTFILES:=Mshichmanhodgesmosfet.oct

all: $(OCTFILES)

%.oct:  %.cc
	mkoctfile $<

clean:
	-rm -f *.o core octave-core *.oct *~
