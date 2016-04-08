all:
	g++ main.cpp func.cpp vtk_cpp.cpp -o gs.exe
clean:
	rm *~ *.exe *.o *.vtk
