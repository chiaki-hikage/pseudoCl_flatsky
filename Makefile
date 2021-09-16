FC2 	= ifort
FC2_MPI	= mpif90
FFLAGS 	= -O2
FFLAGS_MPI  = -DMATRIX_SINGLE -O2 -ip -W0 -WB -qopenmp -fpp -DMPI

PROG = sheardist starmask modemat_square_para modemat_inside_para shearpow deconv collect

OBJPRO = sheardist.o starmask.o modemat_square_para.o modemat_inside_para.o shearpow.o deconv.o collect.o

OBJSUB = sub.o fieldinfo.o

####

.SUFFIXES: .f90 .o .mod

.f90.o: 
	$(FC2) $(FFLAGS) -c $< 

.f90.mod:
	$(FC2) $(FFLAGS) -c $< 

###

shearpow : shearpow.o $(OBJSUB)
	$(FC2) $(FFLAGS) -o $@ $^

starmask : starmask.o $(OBJSUB)
	$(FC2) $(FFLAGS) -o $@ $^ 

sheardist : sheardist.o $(OBJSUB)
	$(FC2) $(FFLAGS) -o $@ $^ 

modemat_square_para.o : modemat_square_para.f90 $(OBJSUB)
	$(FC2_MPI) $(FFLAGS_MPI) -c $<

modemat_square_para : modemat_square_para.o $(OBJSUB)
	$(FC2_MPI) $(FFLAGS_MPI) -o $@ $^ -lm

modemat_inside_para.o : modemat_inside_para.f90 $(OBJSUB)
	$(FC2_MPI) $(FFLAGS_MPI) -c $<

modemat_inside_para : modemat_inside_para.o $(OBJSUB)
	$(FC2_MPI) $(FFLAGS_MPI) -o $@ $^ -lm

deconv : deconv.o $(OBJSUB)
	$(FC2) $(FFLAGS) -o $@ $^ 

collect : collect.o $(OBJSUB)
	$(FC2) $(FFLAGS) -o $@ $^ 

iparam : iparam.o
	$(FC2) $(FFLAGS) -o $@ $^ 

###

shearpow.f90: $(OBJSUB)

starmask.f90: $(OBJSUB) maskinfo.o

sheardist.f90: $(OBJSUB)

modemat_square.f90: $(OBJSUB)

modemat_inside_para.f90: $(OBJSUB)

deconv.f90: $(OBJSUB)

collect.f90: $(OBJSUB)

sub.f90: fieldinfo.o

###

default : $(PROG)

all : $(PROG)

setparam : 
	sh iparam.sh

test : 
	sh test.sh

clean:
	-rm -f $(PROG) iparam *.o *mod core

tidy:
	-rm -f $(PROG) 
