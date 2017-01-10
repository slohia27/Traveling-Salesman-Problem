//SURBHI LOHIA
//SL3893
//DATA STRUCTURES HOMEWORK 6
//MAY 3, 2016

import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;

public class Graph {

	// Keep a fast index to nodes in the map
	private Map<Integer, Vertex> vertexNames;

	/**
	 * Construct an empty Graph with a map. The map's key is the name of a vertex
	 * and the map's value is the vertex object.
	 */
	public Graph() {
		vertexNames = new HashMap<Integer, Vertex>();
	}

	/**
	 * Adds a vertex to the graph. Throws IllegalArgumentException if two vertices
	 * with the same name are added.
	 * 
	 * @param v
	 *          (Vertex) vertex to be added to the graph
	 */
	public void addVertex(Vertex v) {
		if (vertexNames.containsKey(v.name))
			throw new IllegalArgumentException("Cannot create new vertex with existing name.");
		vertexNames.put(v.name, v);
	}

	/**
	 * Gets a collection of all the vertices in the graph
	 * 
	 * @return (Collection<Vertex>) collection of all the vertices in the graph
	 */
	public Collection<Vertex> getVertices() {
		return vertexNames.values();
	}

	/**
	 * Gets the vertex object with the given name
	 * 
	 * @param name
	 *          (String) name of the vertex object requested
	 * @return (Vertex) vertex object associated with the name
	 */
	public Vertex getVertex(String name) {
		return vertexNames.get(name);
	}

	/**
	 * Adds a directed edge from vertex u to vertex v
	 * 
	 * @param nameU
	 *          (String) name of vertex u
	 * @param nameV
	 *          (String) name of vertex v
	 * @param cost
	 *          (double) cost of the edge between vertex u and v
	 */
	public void addEdge(int nameU, int nameV, Double cost) {
		if (!vertexNames.containsKey(nameU))
			throw new IllegalArgumentException(nameU + " does not exist. Cannot create edge.");
		if (!vertexNames.containsKey(nameV))
			throw new IllegalArgumentException(nameV + " does not exist. Cannot create edge.");
		Vertex sourceVertex = vertexNames.get(nameU);
		Vertex targetVertex = vertexNames.get(nameV);
		Edge newEdge = new Edge(sourceVertex, targetVertex, cost);
		sourceVertex.addEdge(newEdge);
	}

	/**
	 * Adds an undirected edge between vertex u and vertex v by adding a directed
	 * edge from u to v, then a directed edge from v to u
	 * 
	 * @param name
	 *          (String) name of vertex u
	 * @param name2
	 *          (String) name of vertex v
	 * @param cost
	 *          (double) cost of the edge between vertex u and v
	 */
	public void addUndirectedEdge(int name, int name2, double cost) {
		addEdge(name, name2, cost);
		addEdge(name2, name, cost);
	}


	/**
	 * Computes the euclidean distance between two points as described by their
	 * coordinates
	 * 
	 * @param ux
	 *          (double) x coordinate of point u
	 * @param uy
	 *          (double) y coordinate of point u
	 * @param vx
	 *          (double) x coordinate of point v
	 * @param vy
	 *          (double) y coordinate of point v
	 * @return (double) distance between the two points
	 */
	public double computeEuclideanDistance(double ux, double uy, double vx, double vy) {
		return Math.sqrt(Math.pow(ux - vx, 2) + Math.pow(uy - vy, 2));
	}

	/**
	 * Computes euclidean distance between two vertices as described by their
	 * coordinates
	 * 
	 * @param u
	 *          (Vertex) vertex u
	 * @param v
	 *          (Vertex) vertex v
	 * @return (double) distance between two vertices
	 */
	public double computeEuclideanDistance(Vertex u, Vertex v) {
		return computeEuclideanDistance(u.x, u.y, v.x, v.y);
	}

	/**
	 * Calculates the euclidean distance for all edges in the map using the
	 * computeEuclideanCost method.
	 */
	public void computeAllEuclideanDistances() {
		for (Vertex u : getVertices())
			for (Edge uv : u.adjacentEdges) {
				Vertex v = uv.target;
				uv.distance = computeEuclideanDistance(u.x, u.y, v.x, v.y);
			}
	}



	// STUDENT CODE STARTS HERE
	
	//words weight/cost/distance more or less interchangeable when used 

	public void generateRandomVertices(int n) {
		vertexNames = new HashMap<Integer, Vertex>(); // reset the vertex hashmap

		//creating random vertex	
		//this will also get us n vertices
		for (int i = 0; i < n; i++) {
			//random allows us get vertices at random coordinates
			int x = (int)((Math.random() * 100) + 1);
			int y = (int)((Math.random() * 100) + 1);
			Vertex vertices = new Vertex(i, x, y);
			addVertex(vertices);
		}

		//creating edges for each vertex
		for (int j = 0; j < n; j++) {
			Vertex vert1 = vertexNames.get(j);
			for (int k = 0; k < n; k++) {
				if (j != k) {
					Vertex vert2 = vertexNames.get(k);
					vertexNames.get(j).adjacentEdges.add(new Edge(vert1, vert2, 0));
				}
			}
		}

		computeAllEuclideanDistances(); // compute distances
	}

	
	public Edge getEdge(int u, int v) {

		Vertex vert1 = vertexNames.get(u);
		//
		for (Edge ver : vert1.adjacentEdges) {
			//going through to find adjacent edges
			//undirected edges
			if (ver.target.name == v)
				return ver;
		}
		return null;
	}

	//need to create a method to progressively add weight of each edge to a total weight
	public int totalWeight(ArrayList<Edge> path) {
		if (path == null)
			return Integer.MAX_VALUE;
		int weight = 0;
		//use for each loop to go through each edge and add each weight
		for (Edge e : path) {
			weight += e.distance;
		}
		return weight;
	}
	
	public List<Edge> nearestNeighborTsp() {
		//this is an arraylist to keep track of the shortest distance from 1 edge to another
		List<Edge> shortestdist = new ArrayList<Edge>();

		//the number of vertices we have in hashmap
		int mapsize = vertexNames.size();
		//creating vertices 
		Vertex vert1 = vertexNames.get((new Random()).nextInt(mapsize));
		Vertex vert2 = vert1;

		//how many vertices have been visited
		int complete = 1;
		boolean visited[] = new boolean[mapsize];
		visited[vert1.name] = true;

		//while we have visited less than the amount of vertices there are
		while (complete < mapsize) {

			Edge currEdge = null;
			for (Edge e : vert2.adjacentEdges) {
				if (!visited[e.target.name] && (currEdge == null
						|| currEdge.distance > e.distance)) {
					currEdge = e;
				}
			}
			shortestdist.add(currEdge);
			vert2 = currEdge.target;
			visited[vert2.name] = true;
			complete++;
		} //undirected edges
		for (Edge e : vert1.adjacentEdges) {
			if (e.target.name == vert2.name) {
				shortestdist.add(e);
				break;
			}
		}

		return shortestdist;
	}
	
	//this swap method is needed for permutations in brute force
	public void swap(int a, int b, int[] permutation) {
		int temp = permutation[a];
		permutation[a] = permutation[b];
		permutation[b] = temp;
	}
	
	//keeping track of weight of edges via arraylist
	public ArrayList<Edge> minweight(ArrayList<Edge> list1, ArrayList<Edge> list2) {
		int totalWeightList1 = totalWeight(list1);
		int totalWeightList2 = totalWeight(list2);
		return (totalWeightList1 < totalWeightList2 ? list1 : list2);
	}
	

	//begin route with 1st city
	public ArrayList<Edge> route(int[] permutation, int begin, int size) {
		if (begin == size - 1) {
			ArrayList<Edge> path = new ArrayList<Edge>();
			for (int i = 0; i < size; i++) {
				int u = permutation[i];
				int v = permutation[(i + 1) % size];
				path.add(getEdge(u, v));
			}
			return path;
		}
		ArrayList<Edge> path = null; 
		for(int i = begin + 1; i< size; i++){
			swap(begin, i, permutation);
			//trying to find least weight/cost/distance path
			path = minweight(path, route(permutation, begin+1, size));
			swap(begin, i, permutation);
		}
		return path;
	}
	
	public List<Edge> bruteForceTsp() {

		int size = vertexNames.size();

		int[] permutation = new int[size];
		for(int i = 0; i< size; i++){
			permutation[i] = i;
		}
		//recursively use this method on each route
		return route(permutation, 0, size);
	}


	// STUDENT CODE ENDS HERE



	/**
	 * Prints out the adjacency list of the graph for debugging
	 */
	public void printAdjacencyList() {
		for (int u : vertexNames.keySet()) {
			StringBuilder sb = new StringBuilder();
			sb.append(u);
			sb.append(" -> [ ");
			for (Edge e : vertexNames.get(u).adjacentEdges) {
				sb.append(e.target.name);
				sb.append("(");
				sb.append(e.distance);
				sb.append(") ");
			}
			sb.append("]");
			System.out.println(sb.toString());
		}
	}
}
