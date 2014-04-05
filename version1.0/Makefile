FC = gfortran
XPSDIR=$(HOME)/XPS
FFLAGS = -Wall -Wextra -march=native -O3 -ffast-math
FFLAGS += $(shell pkg-config --cflags --libs plplotd-f95)

LDFLAGS = -march=native -O3 -ffast-math

LIBS += $(shell pkg-config --cflags --libs plplotd-f95)

PROGRAMS = md

OBJS =
OBJS += plot.o
OBJS += randomseed.o
OBJS += md.o



myprog: $(OBJS)
	$(FC) $(OBJS) -o md $(LIBS) $(LDFLAGS)



%: %.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f95
	$(FC) $(FFLAGS) -c $<
%.o: %.F95
	$(FC) $(FFLAGS) -c $<
%.o: %.f90
	$(FC) $(FFLAGS) -c $<
%.o: %.F90 
	$(FC) $(FFLAGS) -c $<
%.o: %.f
	$(FC) $(FFLAGS) -c $<
%.o: %.F 
	$(FC) $(FFLAGS) -c $<

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD $(PROGRAMS)

veryclean: clean
	rm -f *~ $(PROGRAMS) 
