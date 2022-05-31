#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <vector>
using namespace std;
#define MAXIMUM 9999999;
ifstream in_file("input.txt");
ofstream out_file("file.txt");
ofstream out_kruskal("kruskal.txt");
ofstream out_prim("prim.txt");
class Link {
public:
	int source, destination;
	int weight;
};
class Graph {
public:

	int Vertices_num, Links_num;
	Link* links;

};
Graph* make_graph(int V, int L)
{
	Graph* graph = new Graph;
	graph->Vertices_num = V;
	graph->Links_num = L;

	graph->links = new Link[L];

	return graph;
}
class subset   //a subtree of the graph
{
public:
	int parent;   //arbitrary parent of the subset 
	int rank; //reflects the maximum depth of the tree produced so far
};


//utility function to return the index of the unvisited node with the minimum key value.

int getMinKey(int keys[], bool Included_v[], int Nodes)
{
	int min = MAXIMUM;
	int index;
	for (int i = 0; i < Nodes; i++)	//looping through all vertices.
		if (Included_v[i] == false && keys[i] < min)//it selcets the unvisited node with the least key value.
		{
			min = keys[i], index = i;
		}
	return index;
}

// A utility function to print the 
// constructed MST stored in MST[] 
void printMST(int MST[], int** graph, int V)
{
	int sum = 0;
	cout << "The Minimum Spanning Tree is : \n";
	out_prim << "The Minimum Spanning Tree is : \n";
	cout << "Edge \tWeight\n";
	out_prim << "Edge \tWeight\n";
	for (int i = 1; i < V; i++)
	{
		sum += graph[i][MST[i]]; //summing the costs to get the total minimum cost of the tree.
		cout << MST[i] << " - " << i << " \t" << graph[i][MST[i]] << " \n";		//displaying the parent of the node - the node - the weight of the edge between them.
		out_prim << MST[i] << " - " << i << " \t" << graph[i][MST[i]] << " \n";
	}
	cout << "------------------------------------\n";
	cout << "The minimum total cost = " << sum << endl;
	out_prim << "The minimum total cost = " << sum;
}

//Prim function. It receives the adjacency matrix of the graph and the number of nodes
void primMST(int** graph, int Nodes)
{
	int* MST = new int[Nodes]; //array to store the prim's tree.
	int* keys = new int[Nodes]; //array to store the key values for each node (vertix).
	bool* Included_v = new bool[Nodes]; //array to keep track of the included nodes in the tree by setting its values to true.

	for (int i = 0; i < Nodes; i++) //initializing the key values of all nodes to inifinity and mark them all as unvisited.
	{
		keys[i] = MAXIMUM; Included_v[i] = false;
	}

	keys[0] = 0; //starting by the first node so we make its key value=0 to be started with and mark it as the root of the tree by making its parent node -1.
	MST[0] = -1;
	for (int count = 0; count < Nodes - 1; count++) // looping until all vertices are included in the tree.
	{
		int u = getMinKey(keys, Included_v, Nodes);	//getting the index of the nodes with minimum key value.
		Included_v[u] = true;						//marking it as visited
		for (int v = 0; v < Nodes; v++)				//updating the key values for the adjacency nodes of the node u.
			if (graph[u][v] && Included_v[v] == false && graph[u][v] < keys[v])
				MST[v] = u, keys[v] = graph[u][v];
	}
	printMST(MST, graph, Nodes);					//priniting the tree.
}
void print(int** A, int N)
{
	cout << "using prim\n";
	out_prim << "Results for sample graph in file (file.txt)\n";
	out_prim << "Adjacency Matrix\n";
	cout << "  ";
	out_prim << "  ";
	for (int i = 0; i < N; i++)
	{
		cout << setw(4) << i;
		out_prim << setw(4) << i;
	}
	cout << endl;
	out_prim << endl;
	out_prim << "________________________________________________\n";
	out_prim << "________________________________________________\n";
	for (int r = 0; r < N; r++)
	{
		cout << r << "|";
		out_prim << r << "|";
		for (int c = 0; c < N; c++)
		{
			cout << setw(4) << A[r][c];
			out_prim << setw(4) << A[r][c];
		}cout << endl;
		out_prim << endl;
	}
	out_prim << endl;
	out_prim << "Number of vertices = " << N << endl;
}
// (uses path compression technique) flattenning the tree of subsets by putting smaller depth trees right under the big root
int find(subset subsets[], int x)    // O(log n) by union by rank instead of the simple approach that takes O(n)
{
	// find root and make root as parent of x (path compression)
	int a = x;
	if (subsets[x].parent != a)  // x isn't a root
		subsets[x].parent = find(subsets, subsets[x].parent);

	return subsets[x].parent;
}
//storing ranks is more efficient than storing heights. The height of a node can change 
//during a Find operation, so storing ranks avoids the extra effort of keeping the height correct
void Union(subset subsets[], int a, int b)  // a and b are parents
{
	int a_root = find(subsets, a);
	int b_root = find(subsets, b);
	if (subsets[a_root].rank > subsets[b_root].rank)
		subsets[b_root].parent = a_root;
	else if (subsets[a_root].rank < subsets[b_root].rank)
		subsets[a_root].parent = b_root;

	else
	{
		subsets[b_root].parent = a_root;
		subsets[a_root].rank++;
	}
}

int weight_comp(const void* L1, const void* L2)
{
	Link* x = (Link*)L1;
	Link* y = (Link*)L2;
	return (x->weight > y->weight);
}

void KruskalMST(Graph* graph)
{
	int V = graph->Vertices_num;
	Link* result = new Link[V];

	int r = 0;   //for results
	int se = 0;   //for the increasingly sorted edges 

//  qsort(pointer to base, num of items, sizeof the first pointer, int (*compar)(const void *, const void*)) sorts an array.
	qsort(graph->links, graph->Links_num, sizeof(graph->links[0]), weight_comp); // using quick sort


	subset* subsets = new subset[V];


	for (int v = 0; v < V; ++v)
	{
		subsets[v].parent = v;
		subsets[v].rank = 0;
	}

	while (r < V - 1 && se < graph->Links_num)  // to take only (V - 1) links 
	{

		Link taken_link = graph->links[se];   // taking the link with smallest weight 
		se++;

		int a = find(subsets, taken_link.source);
		int b = find(subsets, taken_link.destination);

		if (a != b)  // no cycle with the MST
		{
			result[r++] = taken_link;
			Union(subsets, a, b);
		}
	}


	out_kruskal << "Number of vertices = " << graph->Vertices_num << endl;
	cout << "Number of non_zero edges = " << graph->Links_num << endl;
	out_kruskal << "Number of non_zero edges = " << graph->Links_num << endl;
	out_kruskal << "-------------------------------------\n";
	out_kruskal << "-------------------------------------\n";

	for (int i = 0; i < graph->Links_num; i++)
	{
		cout << graph->links[i].source << " <-> " << graph->links[i].destination << " weight: " << graph->links[i].weight << endl;
		out_kruskal << graph->links[i].source << " <-> " << graph->links[i].destination << " weight: " << graph->links[i].weight << endl;
	}
	cout << "\n-------------------------------------\n";
	out_kruskal << "\n-------------------------------------\n";

	cout << "The Minimum Spanning Tree is : \n";
	out_kruskal << "The Minimum Spanning Tree is : \n";
	cout << "Edge \tWeight\n";
	out_kruskal << "Edge \tWeight\n";
	int cost = 0;
	for (int i = 0; i < r; i++)
	{
		cost += result[i].weight;
		cout << result[i].source << " - " << result[i].destination << " \t" << result[i].weight << endl;
		out_kruskal << result[i].source << " - " << result[i].destination << " \t" << result[i].weight << endl;
	}
	cout << "------------------------------------\n";
	cout << "The minimum total cost = " << cost;
	out_kruskal << "The minimum total cost = " << cost;
}



void main()
{
	vector <int> from;
	vector <int> to;
	vector <int> weight;
	int V, Wmin, Wmax, L = 0, choice;

	cout << "If you wanna use a random adjacency matrix enter 1,"
		<< " if you want to read it from a file enter any other number\nYour choice:";
	cin >> choice;

	if (choice == 1)
	{
		cout << "Enter number of vertices: ";
		cin >> V;
		cout << "Enter the min and max weights respictivly: ";
		cin >> Wmin >> Wmax;
	}
	else
	{
		in_file >> V;
		cout << "Number of vertices is: " << V << endl;
		cout << "Enter the min and max weights respictivly: ";
		cin >> Wmin >> Wmax;
	}
	int** adj_mat = new int* [V];
	for (size_t i = 0; i < V; i++)
		*(adj_mat + i) = new int[V];
	srand(time(0));
	if (choice == 1)
	{

		for (size_t i = 0; i < V - 1; i++)
			for (size_t j = i + 1; j < V; j++)
			{
				adj_mat[i][j] = rand() % Wmax;
				if (adj_mat[i][j] < Wmin)
					adj_mat[i][j] = 0;
				adj_mat[j][i] = adj_mat[i][j];
				if (adj_mat[i][j] != 0)
					L++;
				from.push_back(i);
				to.push_back(j);
				weight.push_back(adj_mat[i][j]);
			}
		for (size_t i = 0; i < V; i++)
			adj_mat[i][i] = 0;
	}
	else
	{
		for (size_t i = 0; i < V; i++)
			for (size_t j = 0; j < V; j++)
			{
				in_file >> adj_mat[i][j];
				if (adj_mat[i][j] < Wmin)
					adj_mat[i][j] = 0;
				if (adj_mat[i][j] != 0)
					L++;
				from.push_back(i);
				to.push_back(j);
				weight.push_back(adj_mat[i][j]);
			}
		in_file.close();
		L = L / 2;
	}

	print(adj_mat, V);
	cout << "Number of non-zero edges = " << L << endl;
	out_prim << "Number of non-zero edges = " << L << endl;
	cout << "-------------------------------------\n";
	out_prim << "-------------------------------------\n";
	if (choice == 1)
	{
		for (int i = 0; i < from.size(); i++)
		{

			if (weight[i] == 0)
			{

			}
			else
			{
				out_prim << from[i] << " <-> " << to[i] << " weight: " << weight[i] << endl;
				cout << from[i] << " <-> " << to[i] << " weight: " << weight[i] << endl;
			}
		}
	}
	else
	{
		for (int i = 0; i < from.size(); i++)
		{

			if (weight[i] == 0)
			{

			}
			else
			{
				out_prim << from[i] << " <-> " << to[i] << " weight: " << weight[i] << endl;
				cout << from[i] << " <-> " << to[i] << " weight: " << weight[i] << endl;
			}
		}
	}

	cout << endl;
	out_prim << endl;
	cout << "-------------------------------------\n";
	out_prim << "-------------------------------------\n";

	primMST(adj_mat, V);

	cout << "\n\nusing kruskal\n";

	Graph* graph = make_graph(V, L);
	for (size_t i = 0, k = 0; i < V - 1; i++)
		for (size_t j = i + 1; j < V; j++)
		{
			if (adj_mat[i][j] != 0)
			{
				graph->links[k].source = i;
				graph->links[k].destination = j;
				graph->links[k].weight = adj_mat[i][j];
				k++;
			}
		}

	out_file << V << " ";
	for (size_t i = 0; i < V; i++)
		for (size_t j = 0; j < V; j++)
			out_file << adj_mat[i][j] << " ";
	out_file.close();

	out_kruskal << "Results for sample graph in file (file.txt)\n";
	out_kruskal << "Adjacency Matrix\n";
	cout << "  ";
	out_kruskal << "  ";
	for (int i = 0; i < graph->Vertices_num; i++)
	{
		cout << setw(4) << i;
		out_kruskal << setw(4) << i;
	}
	cout << "\n_________________________________________\n";
	out_kruskal << "\n_________________________________________\n";
	for (size_t i = 0; i < graph->Vertices_num; i++)
	{
		cout << i << "|";
		out_kruskal << i << "|";
		for (size_t j = 0; j < graph->Vertices_num; j++)
		{
			cout << setw(4) << adj_mat[i][j];
			out_kruskal << setw(4) << adj_mat[i][j];
		}
		cout << endl;
		out_kruskal << endl;
	}
	out_kruskal << endl;
	KruskalMST(graph);
	out_kruskal.close();
}
