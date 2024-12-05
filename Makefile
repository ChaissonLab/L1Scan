all: l1scan

alnPairs: AlnPairs.cpp
	g++ -O2 AlnPairs.cpp -IWFA/bindings -IWFA -o alnPairs -L WFA/lib -lwfacpp

l1scan: L1Scan.cpp
	g++ -O2 L1Scan.cpp -o l1scan -Iseqan/include  -I$(CONDA_PREFIX)/include -L$(CONDA_PREFIX)/lib -lhts 


align2: Align2.cpp
	g++ -std=c++17 -Wno-deprecated-declarations Align2.cpp -Iseqan/include -L$(CONDA_PREFIX)/lib -o $@
