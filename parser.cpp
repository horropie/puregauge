#include <iostream>
#include <typeinfo>

void setlattice(int (&dimv)[3]){

	/** We have to check if the input is of integer type. Note that this (unfortunately) accepts things like "10abc", as its accepts the number and ignores the rest. 
    One could split a header isinteger.h from this to make it simpler.**/
	int isinteger;

	std::cout << "What is the size of the lattice in x-direction?" << std::endl;
    std::cin >> isinteger;
    while(std::cin.fail()) {
        std::cout << "Please enter an integer value." << std::endl;
        std::cin.clear();
        std::cin.ignore(256,'\n');
        std::cin >> isinteger;
    }
    isinteger >> dimv[0];

	std::cout << "What is the size of the lattice in y-direction?" << std::endl;
    std::cin >> isinteger;
    while(std::cin.fail()) {
        std::cout << "Please enter an integer value." << std::endl;
        std::cin.clear();
        std::cin.ignore(256,'\n');
        std::cin >> isinteger;
    }
	isinteger >> dimv[1];

	std::cout << "What is the size of the lattice in z-direction?" << std::endl;
    std::cin >> isinteger;
    while(std::cin.fail()) {
        std::cout << "Please enter an integer value." << std::endl;
        std::cin.clear();
        std::cin.ignore(256,'\n');
        std::cin >> isinteger;
    }
	isinteger >> dimv[2];

}

/**
int main() {
	int dimv[3];

	setlattice(dimv);	

	return 0;
}
**/
