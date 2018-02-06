all:
	gfortran -c get_all_fields_mod.f90
	gfortran -c halo_exchange_mod.f90
	gfortran -c sea_breeze_diag.f90
	gfortran -c dummy_model.f90
	gfortran get_all_fields_mod.o halo_exchange_mod.o sea_breeze_diag.o dummy_model.o -o dummy_model
clean:
	rm -f *.o *.mod dummy_model
