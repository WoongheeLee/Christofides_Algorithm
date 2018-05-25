import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeMap;

import org.jgrapht.alg.cycle.HierholzerEulerianCycle;
import org.jgrapht.alg.spanning.KruskalMinimumSpanningTree;
import org.jgrapht.graph.AsSubgraph;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleWeightedGraph;
import org.jgrapht.graph.WeightedMultigraph;
import org.jgrapht.traverse.DepthFirstIterator;

public class Christofides {
	public static ArrayList<Integer> get2ApproximationTSP(SimpleWeightedGraph<Integer,DefaultWeightedEdge> G) {
		// STEP 1: Get Minimum Spanning Tree
		SimpleWeightedGraph<Integer,DefaultWeightedEdge> minimum_spanning_tree = getMinimumSpanningTree(G);
				
		ArrayList<Integer> path = new ArrayList<Integer>();
		HashSet<Integer> temp_set = new HashSet<Integer>();

		// STEP 2: DFS and Skip a Node if it is visited.
		for (DepthFirstIterator<Integer,DefaultWeightedEdge> i = new DepthFirstIterator<>(minimum_spanning_tree); i.hasNext(); ) {
			int v = i.next();
			if (!temp_set.contains(v)) {
				temp_set.add(v);
				path.add(v);
			}
		}
		
		return path;
	}

	// STEP 1: Get Minimum Spanning Tree
	private static SimpleWeightedGraph<Integer,DefaultWeightedEdge> getMinimumSpanningTree(SimpleWeightedGraph<Integer, DefaultWeightedEdge> G) {
		SimpleWeightedGraph<Integer, DefaultWeightedEdge> minimum_spanning_tree = new SimpleWeightedGraph<Integer, DefaultWeightedEdge>(DefaultWeightedEdge.class);
		
		for (Iterator<Integer> i = G.vertexSet().iterator(); i.hasNext(); )
			minimum_spanning_tree.addVertex(i.next());
		
		KruskalMinimumSpanningTree<Integer, DefaultWeightedEdge> kruskal = new KruskalMinimumSpanningTree<Integer, DefaultWeightedEdge>(G);
		for (Iterator<DefaultWeightedEdge> i = kruskal.getSpanningTree().iterator(); i.hasNext(); ) {
			DefaultWeightedEdge e = i.next();
			int v1 = G.getEdgeSource(e);
			int v2 = G.getEdgeTarget(e);
			double w = G.getEdgeWeight(e);
			
			minimum_spanning_tree.addEdge(v1, v2);
			DefaultWeightedEdge temp_e = minimum_spanning_tree.getEdge(v1, v2);
			minimum_spanning_tree.setEdgeWeight(temp_e, w);
		}
		
		return minimum_spanning_tree;
	}
	
	// STEP 2-1: Get Odd Degree Set
	private static HashSet<Integer> getOddDegreeVertices(
			SimpleWeightedGraph<Integer, DefaultWeightedEdge> minimum_spanning_tree) {

		HashSet<Integer> odd_degree_set = new HashSet<Integer>();
		
		for (Iterator<Integer> i = minimum_spanning_tree.vertexSet().iterator(); i.hasNext(); ) {
			int v = i.next();
			if (minimum_spanning_tree.degreeOf(v)%2 != 0)
				odd_degree_set.add(v);
		}
		
		return odd_degree_set;
	}
	
	// STEP 3: Get Minimum Weight Perfect Mathing from the Subgraph
	private static SimpleWeightedGraph<Integer,DefaultWeightedEdge> getMinimumWeightedPerfectMatching(AsSubgraph<Integer, DefaultWeightedEdge> subgraph_G) {
		SimpleWeightedGraph<Integer,DefaultWeightedEdge> min_weight_perfect_matching = new SimpleWeightedGraph<Integer,DefaultWeightedEdge>(DefaultWeightedEdge.class);
		TreeMap<Double, ArrayList<DefaultWeightedEdge>> edge_map = new TreeMap<Double,ArrayList<DefaultWeightedEdge>>();
		
		// Sort Edges by Weight
		for (Iterator<DefaultWeightedEdge> i = subgraph_G.edgeSet().iterator(); i.hasNext(); ) {
			DefaultWeightedEdge e = i.next();
			double w = subgraph_G.getEdgeWeight(e);
			ArrayList<DefaultWeightedEdge> temp_e_list = null;
			
			if (edge_map.containsKey(w)) 
				temp_e_list = edge_map.get(w);
			else
				temp_e_list = new ArrayList<DefaultWeightedEdge>();
			temp_e_list.add(e);
			
			edge_map.put(w, temp_e_list);
		}
		
		// Get Minimum Weight Perfect Matching
		for (Iterator<Double> i = edge_map.keySet().iterator(); i.hasNext(); ) {
			double w = i.next();
			ArrayList<DefaultWeightedEdge> e_list = edge_map.get(w);
			
			for (DefaultWeightedEdge e : e_list) {
				int v1 = subgraph_G.getEdgeSource(e);
				int v2 = subgraph_G.getEdgeTarget(e);
				
				if (min_weight_perfect_matching.containsVertex(v1) || 
						min_weight_perfect_matching.containsVertex(v2))
					continue;
				else {
					min_weight_perfect_matching.addVertex(v1);
					min_weight_perfect_matching.addVertex(v2);
					min_weight_perfect_matching.addEdge(v1, v2);
					DefaultWeightedEdge temp_e = min_weight_perfect_matching.getEdge(v1, v2);
					min_weight_perfect_matching.setEdgeWeight(temp_e, w);
				}
			}
			
			if (min_weight_perfect_matching.vertexSet().size() == 
					subgraph_G.vertexSet().size()) break;
		}
		
		return min_weight_perfect_matching;
	}
	
	// STEP 4: Combine Perfect Matching and Minimum Spanning Tree
	private static WeightedMultigraph<Integer,DefaultWeightedEdge> getMultiGraph(
			SimpleWeightedGraph<Integer,DefaultWeightedEdge> minimum_spanning_tree, 
			SimpleWeightedGraph<Integer,DefaultWeightedEdge> min_weight_perfect_matching) {
		
		WeightedMultigraph<Integer,DefaultWeightedEdge> multi_graph =
				new WeightedMultigraph<Integer,DefaultWeightedEdge>(DefaultWeightedEdge.class);
		
		// Add Vertex
		for (Iterator<Integer> i = minimum_spanning_tree.vertexSet().iterator(); i.hasNext(); )
			multi_graph.addVertex(i.next());
		for (Iterator<Integer> i = min_weight_perfect_matching.vertexSet().iterator(); i.hasNext(); ) {
			int v = i.next();
			if (!multi_graph.containsVertex(v))
				multi_graph.addVertex(v);
		}
		
		// Add Edge
		for (Iterator<DefaultWeightedEdge> i = minimum_spanning_tree.edgeSet().iterator(); i.hasNext(); ) {
			DefaultWeightedEdge e = i.next();
			int v1 = minimum_spanning_tree.getEdgeSource(e);
			int v2 = minimum_spanning_tree.getEdgeTarget(e);
			double w = minimum_spanning_tree.getEdgeWeight(e);
			
			multi_graph.addEdge(v1, v2);
			DefaultWeightedEdge temp_e = multi_graph.getEdge(v1, v2);
			multi_graph.setEdgeWeight(temp_e, w);
		}
		for (Iterator<DefaultWeightedEdge> i = min_weight_perfect_matching.edgeSet().iterator(); i.hasNext(); ) {
			DefaultWeightedEdge e = i.next();
			int v1 = min_weight_perfect_matching.getEdgeSource(e);
			int v2 = min_weight_perfect_matching.getEdgeTarget(e);
			double w = min_weight_perfect_matching.getEdgeWeight(e);
			
			multi_graph.addEdge(v1, v2);
			DefaultWeightedEdge temp_e = multi_graph.getEdge(v1, v2);
			multi_graph.setEdgeWeight(temp_e, w);
		}

		return multi_graph;
	}
	
	// STEP 5: Get Eulerian Circuit
	private static ArrayList<Integer>
		getEulerTour(WeightedMultigraph<Integer,DefaultWeightedEdge> multi_graph) {
		
		HierholzerEulerianCycle<Integer,DefaultWeightedEdge> h = new HierholzerEulerianCycle<Integer,DefaultWeightedEdge>();
		
		ArrayList<Integer> euler_tour = new ArrayList<Integer>();
		
		int v = h.getEulerianCycle(multi_graph).getStartVertex();
		
		for (DefaultWeightedEdge e : h.getEulerianCycle(multi_graph).getEdgeList()) {
			int s = multi_graph.getEdgeSource(e);
			int t = multi_graph.getEdgeTarget(e);
			if (v == s) {
				euler_tour.add(s);
				v = t;
			} else {
				euler_tour.add(t);
				v = s;
			}
		}
		euler_tour.add(v);
		
		return euler_tour;
	}
	
	// STEP 6: Remove Repeated Vertices (Finding Hamilton Path)
	private static ArrayList<Integer> getHamiltonianPath(ArrayList<Integer> euler_tour) {
		// This method fits only in the Christofides class.
		ArrayList<Integer> hamiltonian = new ArrayList<Integer>();
		HashSet<Integer> temp_set = new HashSet<Integer>();
		
		for (int t : euler_tour) {
			if (!temp_set.contains(t)) {
				temp_set.add(t);
				hamiltonian.add(t);
			}
		}
		hamiltonian.add(euler_tour.get(0));
		
		return hamiltonian;
	}
	
	public static ArrayList<Integer> getTSP(SimpleWeightedGraph<Integer,DefaultWeightedEdge> G) {
		ArrayList<Integer> path = new ArrayList<Integer>();

		// Christofides' Algorithm 
		
		// STEP 1: Get Minimum Spanning Tree
		SimpleWeightedGraph<Integer,DefaultWeightedEdge> minimum_spanning_tree = getMinimumSpanningTree(G);
		
		// STEP 2-1: Get Odd Degree Set
		HashSet<Integer> odd_degree_set = getOddDegreeVertices(minimum_spanning_tree);
		
		// STEP 2-2: Inducing Subgraph with Odd Degree Set
		AsSubgraph<Integer,DefaultWeightedEdge> subgraph_G = new AsSubgraph(G, odd_degree_set, G.edgeSet());
		
		// STEP 3: Get Minimum Weight Perfect Mathing from the Subgraph
		SimpleWeightedGraph<Integer,DefaultWeightedEdge> min_weight_perfect_matching =
				getMinimumWeightedPerfectMatching(subgraph_G);
		
		// STEP 4: Combine Perfect Matching and Minimum Spanning Tree
		WeightedMultigraph<Integer,DefaultWeightedEdge> multi_graph = 
				getMultiGraph(minimum_spanning_tree, min_weight_perfect_matching);
		
		// STEP 5: Get Eulerian Circuit
		ArrayList<Integer> euler_tour = getEulerTour(multi_graph);
		
		// STEP 6: Remove Repeated Vertices (Finding Hamilton Path)
		path = getHamiltonianPath(euler_tour);
		
		return path;
	}
		
}
