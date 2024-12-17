all: l1scan

alnPairs: AlnPairs.cpp
	g++ -O2 AlnPairs.cpp -IWFA/bindings -IWFA -o $@ -L WFA/lib -lwfacpp

l1scan: L1Scan.cpp
	g++ -DSEQAN_ENABLE_TESTING -O2 L1Scan.cpp -o $@ -Iseqan/include  -I$(CONDA_PREFIX)/include -L$(CONDA_PREFIX)/lib -lhts 


minimalTest: MinimalTest.cpp
	g++ -DSEQAN_ENABLE_TESTING -g MinimalTest.cpp -o $@ -Iseqan/include 


align2: Align2.cpp
	g++ -std=c++17 -Wno-deprecated-declarations Align2.cpp -Iseqan/include -L$(CONDA_PREFIX)/lib -o $@
