
#include <iostream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <vector>

#include "smoothsort.h"

void test_stl_sorting() {
	srand(time(NULL));
	std::vector<double> to_sort;
	to_sort.reserve(10);
	for(unsigned int i=0; i<10; i++) {
		to_sort.push_back((int)(100*(double)rand()/(double)RAND_MAX));
	}
	
	sorted_ptr_arr<double> sorted_data(to_sort);
	int last = sorted_data[0];
	bool not_sorted = false;
	for(sorted_ptr_arr<double>::iterator it=sorted_data.begin(); it != sorted_data.end(); ++it) {
		if(last > *it) { not_sorted = true; }
		std::cout << *it << "\t";
	}
	std::cout << std::endl;
	if(not_sorted) { std::cout << "# NOT SORTED" << std::endl; }
}

void test_smoothsort(unsigned int N_lists, unsigned int list_max) {
	srand(time(NULL));
	
	std::cout << "Sorting " << N_lists << " lists of at most " << list_max << " elements each..."; 
	for(unsigned int k=0; k<N_lists; k++) {
		unsigned int N = list_max;//= (int)((double)list_max*(double)rand()/(double)RAND_MAX)+1;
		int* test_input = new int[N];
		int* input_copy = new int[N];
		for(unsigned int i=0; i<N; i++) { test_input[i] = (int)(100*(double)rand()/(double)RAND_MAX); input_copy[i] = test_input[i]; }
		
		//std::cout << std::endl;
		//for(unsigned int i=0; i<N; i++) { std::cout << test_input[i] << "\t"; }
		//std::cout << std::endl << std::endl;
		
		sorted_ptr_arr<int> sorted_data(&test_input[0], N);
		int last = sorted_data[0];
		bool not_sorted = false;
		for(sorted_ptr_arr<int>::iterator it=sorted_data.begin(); it != sorted_data.end(); ++it) {
			if(last > *it) { not_sorted = true; }
			//std::cout << *it << "\t";
		}
		//std::cout << std::endl;
		if(not_sorted) { std::cout << "# NOT SORTED" << std::endl; }
		//std::cout << std::endl;
		
		
		delete[] test_input;
		delete[] input_copy;
	}
	std::cout << "Done." << std::endl;
}

void test_resorting(unsigned int N_trials, unsigned int list_size) {
	srand(time(NULL));
	
	int* test_input = new int[list_size];
	for(unsigned int i=0; i<list_size; i++) { test_input[i] = (int)(100*(double)rand()/(double)RAND_MAX); }
	sorted_ptr_arr<int> sorted_data(&test_input[0], list_size);
	
	std::cout << "Sorting nearly-ordered list of " << list_size << " integers " << N_trials << " times..."; 
	for(unsigned int k=0; k<N_trials; k++) {
		for(unsigned int i=0; i<list_size; i++) { test_input[i] += (int)(3*(double)rand()/(double)RAND_MAX); }
		sorted_data.sort();
	}
	std::cout << "Done." << std::endl;
	
	delete[] test_input;
}

void print_usage(char *argv_0) {
	std::cout << "Usage: " << argv_0 << " n_1 n_2 ... n_N\tSort {n_i}." << std::endl;
	std::cout << "Usage: " << argv_0 << " --stl_test\tTest stl-container sorting." << std::endl;
	std::cout << "Usage: " << argv_0 << " --test N max_size\tSort N lists of size less than max_size." << std::endl;
	std::cout << "Usage: " << argv_0 << " --resort N size\tSort nearly-ordered list of given size N times." << std::endl;
}

int main(int argc, char **argv) {
	if(argc > 1) {
		if(std::string(argv[1]) == "--test") {
			if(argc == 4) {
				unsigned int N_lists = (unsigned int)atoi(argv[2]);
				unsigned int list_max = (unsigned int)atoi(argv[3]);
				test_smoothsort(N_lists, list_max);
			} else {
				print_usage(argv[0]);
			}
		} else if(std::string(argv[1]) == "--resort") {
			if(argc == 4) {
				unsigned int N_lists = (unsigned int)atoi(argv[2]);
				unsigned int list_max = (unsigned int)atoi(argv[3]);
				test_resorting(N_lists, list_max);
			} else {
				print_usage(argv[0]);
			}
		} else if(std::string(argv[1]) == "--stl_test") {
			if(argc == 2) {
				test_stl_sorting();
			} else {
				print_usage(argv[0]);
			}
		} else {
			double* test_input = new double[argc-1];
			for(int i=0; i<argc-1; i++) { test_input[i] = atof(argv[i+1]); }
			smoothsort(test_input, (size_t)(argc-1));
			for(int i=0; i<argc-1; i++) { std::cout << (i == 0 ? "" : "\t") << test_input[i]; }
			std::cout << std::endl;
		}
	} else {
		print_usage(argv[0]);
		return 0;
	}
	
	return 0;
}
