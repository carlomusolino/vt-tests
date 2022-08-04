CXX=g++
CXX_FLAGS=-std=c++14 -g -Wall -I${I_MPI_ROOT}/intel64/include -I${VT_ROOT}/include -I${CHKPT_ROOT}/include -I${DETECTOR_ROOT}/include
LD_FLAGS=-L${I_MPI_ROOT}/intel64/lib -L${VT_ROOT}/lib -L${CHKPT_ROOT}/lib -rdynamic -lEngFormat-Cpp -lfort -lbrotlicommon -lbrotlidec -lbrotlienc -lmpicxx -lcheckpoint -lvt-debug -lfmt-vtd -lz -ldl -lpthread

all: exe

main.o: main.cc include/mesh.hh include/node.hh include/mesh_array.hh include/mesh_data.hh
	$(CXX) $(CXX_FLAGS) -c main.cc -o main.o

exe: main.o
	$(CXX) main.o $(LD_FLAGS) -o main.out 

clean:
	rm main.o main.out