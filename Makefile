
CC=g++ --std=c++11
CFLAGS=-W -Wall -ansi -pedantic
LDFLAGS=
EXEC=compile

all: $(EXEC)
	./main

$(EXEC):demo/main.cpp L2Fsim/flight_zone/flat_thermal_soaring_zone.cpp L2Fsim/flight_zone/thermal/std_thermal.cpp
	$(CC) -o main demo/main.cpp L2Fsim/flight_zone/flat_thermal_soaring_zone.cpp L2Fsim/flight_zone/thermal/std_thermal.cpp -I.

energy:
	python3 python_plot/energy.py

traj:
	python3 python_plot/traj.py

angle:
	python3 python_plot/angle.py

therm:
	python3 python_plot/therm_zslice.py

thermscat:
	python3 python_plot/thermalscatter.py

clean:
	rm -f DATA/data_plane.txt
	rm -f DATA/wind.txt
	rm -f DATA/config.txt
	rm -f DATA/energy.txt
	rm -f main

help:
	@echo “all : run ./main”
	@echo “therm : run therm.py”
	@echo “traj : run traj.py”
	@echo “energy : run energy.py”
	@echo “angle : run angle.py”
	@echo “clean : make clean ...”
