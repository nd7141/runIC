// Calculate Independent Cascade size
// Input: <filename of the weighted graph> <number of vertices> <filename of seeds> <max_iter>
// Example: hep15233.txt 15233 seeds.txt 100

#include <iostream>
#include <fstream> // for reading files
#include <cstdlib> // for atoi, rand
#include <map>
#include <tuple> // for tie
#include <sys/time.h> // for gettimeofday
#include "dirent.h" 				// for reading files from directory
#include <stdio.h> // printf()
#include <stdlib.h> // exit()
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iomanip>
#include <sstream>
#include <string>


#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

using namespace boost;
using namespace std;

typedef property<edge_weight_t, double> Weight;
typedef adjacency_list <vecS, vecS, undirectedS, no_property, Weight > WeightedGraph;
typedef adjacency_list < vecS, vecS, undirectedS> DeterministicGraph;
typedef typename boost::graph_traits<DeterministicGraph>::vertex_descriptor Vertex;
typedef typename graph_traits<WeightedGraph>::edge_descriptor Edge;
typedef typename property_map <WeightedGraph, edge_weight_t>::type myWeight;


typename graph_traits <WeightedGraph>::out_edge_iterator out1, out2;
myWeight weight;


// return execution time
double print_time(struct timeval &start, struct timeval &end){
	double t1=start.tv_sec+(start.tv_usec/1000000.0);
	double t2=end.tv_sec+(end.tv_usec/1000000.0);
	return t2-t1;
}

// read Graph with probabilities from file
void buildList (WeightedGraph& G, string dataset_f) {
	typedef typename WeightedGraph::edge_property_type Weight;
	typename graph_traits<WeightedGraph>::edge_iterator ei, eiend;

	map<int, int> m;

	ifstream infile(dataset_f.c_str());
	if (infile==NULL){
		cout << "Unable to open the input file\n";
	}
	int u, v;
	double p;
	Edge e;
	bool inserted;
	int mapped=0;
	int edge_count=0;

	while (infile >> u >> v >> p){
		tie(e,inserted) = add_edge(u, v, Weight(p), G);
		if (!inserted) {
			cout << "Unable to insert edge\n";
		}
		edge_count++;
	}
//	cout <<"E = "<<edge_count<< " V = "<< mapped <<endl;
}

// read seed nodes from the file
void getSeeds(vector<int>& S, string seeds_f) {
	ifstream infile(seeds_f.c_str());
	if (infile==NULL){
		cout << "Unable to open the input file\n";
	}
	int node;
	while (infile >> node) {
		S.push_back(node);
	}
}

// --------------- MAIN ----------------------//
int main(int argc, char* argv[]) {
	srand(time(NULL));

	struct timeval ex_start, ex_finish;
	gettimeofday(&ex_start, NULL);

	// read parameters from command-line
	const string dataset_f = argv[1]; // filename of the dataset
	const int V = atoi(argv[2]); // number of nodes
	const string seeds_F = argv[3]; // filename of seeds
	const int I = atoi(argv[4]);         // number of MC simulations

	cout << "Graph: " << dataset_f << " I: " << I << endl;

	// read graph from the file
	WeightedGraph G(V);

	struct timeval t_start, t_finish;
	gettimeofday(&t_start, NULL);
	buildList(G, "datasets/" + dataset_f);
	gettimeofday(&t_finish, NULL);
	cout << "* Read original graph: " << print_time(t_start, t_finish) << " sec." << endl;

	// for every file containing seed set per line calculate spread
	string seeds_dir = "seeds/Sparsified_seeds/";
	for (int K=10; K<=100; K += 10) {
		string seeds_f = seeds_dir + seeds_F + lexical_cast<string>(K) + ".txt";
		cout << seeds_f << endl;

		// read seeds from the file
		ifstream infile(seeds_f);
		string line;
		int k;
		while(getline(infile, line)) {
			vector<int> S;
			int count = 0;
			tokenizer<> tok(line);
			gettimeofday(&t_start, NULL);
			for(tokenizer<>::iterator beg=tok.begin(); beg!=tok.end();++beg) {
				int v=lexical_cast<int>(*beg);
				// the first number in line is seed set size k
				if (count == 0) {
					k = v;
					count++;
				}
				// all others are nodes
				else {
					S.push_back(v);
					count++;
				}
			}
			gettimeofday(&t_finish, NULL);
//			cout << "* Read " << S.size() << " seeds: " << print_time(t_start, t_finish) << " sec." << endl;

			// calculate average cascade size
			double cascade_size = 0;
			int iter = 0;
			struct timeval iter_start, iter_stop;
			map<int, bool> activated;
			vector<int> T;
			gettimeofday(&iter_start, NULL);
			while (iter < I) {
				// activate seeds
				for (int i = 0; i < V; ++i) {
					activated[i] = false;
				}
				for (vector<int>::iterator it = S.begin(); it != S.end(); ++it) {
					activated[*it] = true;
					T.push_back(*it);
				}

				// activate new nodes
				vector<int>::size_type ix = 0; // activated node index
				while (ix < T.size()) {
					Vertex u = vertex(T[ix], G);
					int power = 0; // number of nodes u activates
					for (tie(out1, out2) = out_edges(u, G); out1 != out2; ++out1) {
						Vertex v = target(*out1, G);
						if (!activated[v]) {
							double random = (double) rand()/RAND_MAX;
							double p = get(weight, *out1);
							if (random <= p) {
								activated[v] = true;
								T.push_back(v);
								power++;
							}
						}
					}
					ix++;
				}

				cascade_size += T.size();
				iter++;
				gettimeofday(&iter_stop, NULL);
//				cout << iter << ": " << print_time(iter_start, iter_stop) << " sec. Cascade: " << T.size() << endl;
				T.clear();
			}
			gettimeofday(&iter_stop, NULL);
//			cout << "* Average iteration: " << print_time(iter_start, iter_stop)/I << " sec." << endl;
//			cout << "* All iterations: " << print_time(iter_start, iter_stop) << " sec." << endl;

			cascade_size = (double)cascade_size/(double)I;
			cout << "k:" << k << " Average cascade size: " << cascade_size << endl;

			// append cascade size to file
			std::ofstream myfile;
			string output = "spread/K" + lexical_cast<string>(K) + ".txt";
			myfile.open(output.c_str(), ios_base::app);
			myfile << k << " " << cascade_size << endl;
			myfile.close();

			// clear seed set
			S.clear();
		}
	}
	gettimeofday(&ex_finish, NULL);
	cout << "* Execution time: " << print_time(ex_start, ex_finish) << " sec." << endl;
	return 0;
}
