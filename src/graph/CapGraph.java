/**
 *
 */
package graph;
import util.GraphLoader;

import java.util.*;


/**
 * @author Your name here.
 *
 * For the warm up assignment, you must implement your Graph in a class
 * named CapGraph.  Here is the stub file.
 *
 */
public class CapGraph implements Graph {

	/* (non-Javadoc)
	 * @see graph.Graph#addVertex(int)
	 * AdjListsMap is used for week 1 (warm-up) and easy question (get dominating set) for week 3
	 * edges and nodes are used for harder question (detecting sub-communities using Garvin - Newman method)
	 *
	 */
	private HashMap<Integer, HashSet<Integer>> adjListsMap = new HashMap<Integer, HashSet<Integer>>();
	private HashSet<Edge> edges = new HashSet<Edge>();
	private HashMap<Integer, Node> nodes = new HashMap<Integer,Node>();
	private long maxBetweenness = Math.round(Double.POSITIVE_INFINITY);

	@Override
	/* Common for all weeks
	 * (non-Javadoc)
	 * @see graph.Graph#addVertex(int)
	 */
	public void addVertex(int num) {
		// TODO Auto-generated method stub
		if (adjListsMap.get(num) == null){
			HashSet<Integer> neighbors = new HashSet<Integer>();
			adjListsMap.put(num,  neighbors);
			nodes.put(num,new Node(num));
		}

	}

	/* (non-Javadoc)
	 * @see graph.Graph#addEdge(int, int)
	 * Common for all weeks
	 */
	@Override
	public void addEdge(int from, int to) {
		// TODO Auto-generated method stub
		adjListsMap.get(from).add(to);
		Node n1 = nodes.get(from);
		Node n2 = nodes.get(to);
		Edge edge = new Edge(n1,n2);
		edges.add(edge);
		n1.addEdge(edge);
	}

	/* (non-Javadoc)
	 * @see graph.Graph#getEgonet(int)
	 * Week 1 : warm-up assignment
	 */
	@Override
	public Graph getEgonet(int center) {
		// TODO Auto-generated method stub
		CapGraph egoGraph = new CapGraph();
		// Add the immediate neighbors of the center
		HashSet neighbors = adjListsMap.get(center);
		if (neighbors == null){
			return null;
		} else
		{
			egoGraph.addVertex(center);
			Iterator<Integer> iter = neighbors.iterator();
			while (iter.hasNext()){
				int friend = iter.next();
				egoGraph.addEdge(center,friend);
				egoGraph.addVertex(friend);
				// find any edges between neighbors of the center
				HashSet egoFriends = adjListsMap.get(friend);
				if (egoFriends != null ){
					Iterator<Integer> egoIter = egoFriends.iterator();
					while (egoIter.hasNext()){
						int egoFriendNext = egoIter.next();
						if (neighbors.contains(egoFriendNext)){
							egoGraph.addEdge(friend,egoFriendNext);

						}
					}
				}


			}

		}


		return egoGraph;
	}

	/* (non-Javadoc)
	 * @see graph.Graph#getSCCs()
	 * Week 1 : warm-up assignment
	 */
	@Override
	public List<Graph> getSCCs() {
		// TODO Auto-generated method stub
		List<Graph> outG = new ArrayList<Graph>();
		Stack<Integer> vertices = new Stack<Integer>();
		for (int v : adjListsMap.keySet()){
			vertices.push(v);
		}
		Stack<Integer> step1_vertices = dfs(this,vertices);
		CapGraph GT = transposeGraph();
		dfs2(GT,step1_vertices,outG);


		return outG;
	}

	/* Depth first search - First time
	 * Week 1 : warm-up assignment
	 */
	public Stack<Integer> dfs(Graph G, Stack<Integer> vertices){
		HashSet<Integer> visited = new HashSet<Integer>();
		Stack<Integer> finished = new Stack<Integer>();
		while (!vertices.empty()){
			int v = vertices.pop();
			if (!visited.contains(v)){
				dfsVisit(G,v,visited,finished);
			}
		}
		return finished;
	}

	/* Depth first search - 2nd time
	 * Week 1 : warm-up assignment
	 */
	public void dfs2(CapGraph GT, Stack<Integer> vertices, List<Graph> outG){
		HashSet<Integer> visited = new HashSet<Integer>();
		while (!vertices.empty()){
			int v = vertices.pop();
			CapGraph sccGraph = new CapGraph();
			if (!visited.contains(v)){
				dfsVisit2(GT,v,visited,sccGraph);
			}
			if (sccGraph.getNumVertices() > 0){
				outG.add(sccGraph);
			}
		}
	}

	/* dfsvisit first time
	 * Week 1 : warm-up assignment
	 */
	public void dfsVisit(Graph G, int v, HashSet<Integer> visit, Stack<Integer> finish){
		visit.add(v);
		HashSet<Integer> neighbors = getNeighbors(v);
		if (neighbors!=null){
			for (int n : neighbors){
				if (!visit.contains(n)){
					dfsVisit(G,n,visit,finish);
				}
			}
			finish.push(v);
		}
	}

	/*
	 * Week 1 : warm-up assignment
	 * dfsvisit - 2nd time
	 */
	public void dfsVisit2(CapGraph GT, int v, HashSet<Integer> visit, Graph scc){
		visit.add(v);
		HashSet<Integer> neighbors = GT.getNeighbors(v);
		if (neighbors!=null){
			for (int n : neighbors){
				if (!visit.contains(n)){
					dfsVisit2(GT,n,visit,scc);
				}
			}
			scc.addVertex(v);
		}
	}
/*
 * Week 1 : warm-up assignment and week 3 (Capstone project: easy question)
 */
	public HashSet<Integer> getNeighbors(int v){
		HashSet<Integer> neighbors = adjListsMap.get(v);
		return neighbors;
	}

	/*
	 * Capstone project: harder question (Edge betweenness)
	 */

	public HashSet<Node> getNodeList(){
		HashSet<Node> nodeList = new HashSet<Node>();
		for (int id: nodes.keySet()){
			nodeList.add(nodes.get(id));
		}
		return nodeList;
	}

	/*
	 * Capstone project: harder question (Edge betweenness)
	 */

	public HashSet<Edge> getEdgeList(){
		return edges;
	}

	/*
	 * Week 1 - warm-up assignment
	 */

	public CapGraph transposeGraph(){
		CapGraph GT = new CapGraph();
		for (int i : adjListsMap.keySet()){
			HashSet<Integer> neighbors = getNeighbors(i);
			Iterator<Integer> iter = neighbors.iterator();
			while (iter.hasNext()){
				int friend = iter.next();
				GT.addVertex(friend);
				GT.addEdge(friend,i);
			}
			GT.addVertex(i);
		}
		return GT;
	}
	/* (non-Javadoc)
	 * @see graph.Graph#exportGraph()
	 * week 1 - warm-up assigment
	 */
	@Override
	public HashMap<Integer, HashSet<Integer>> exportGraph() {
		// TODO Auto-generated method stub
		return adjListsMap;
	}

	/*
	 * Capstone project: print the whole graph
	 */

	public String printGraph(){
		String s = "";
		s += " (size " + getNumVertices() + "+" + getNumEdges() + " Graph nodes and edges ):";

		for (int k : adjListsMap.keySet()){
			s += "\n\t" + k + " ==> ";
			for (int v : adjListsMap.get(k)){
				s += "\t" + v;
			}
		}
		return s;

	}

/*
 * Capstone project: Get the total number of vertices in the graph
 */
	public int getNumVertices()
	{
		return nodes.keySet().size();
	}

	/*
	 * Capstone project: Get the total number of Edges in the graph
	 */

	public int getNumEdges()
	{
		return edges.size();
	}

	// Project week 3 - assignment easier question: To find the dominating set of nodes in a  network

	public Set<Integer> getDominatingSet() {
		//output set of the dominating set of nodes
		Set<Integer> dominatingSet = new HashSet<Integer>();
		Set<Integer> uncoveredNodes = new HashSet<Integer>();

		//Mark all the nodes as uncovered

		Set<Integer> keyNodes = adjListsMap.keySet();
		Iterator<Integer> nodeIter = keyNodes.iterator();
		while (nodeIter.hasNext()){
			int i = nodeIter.next();
			Set<Integer> uncoveredFriends = adjListsMap.get(i);
			Iterator<Integer> friendIter = uncoveredFriends.iterator();
			while (friendIter.hasNext()){
				uncoveredNodes.add(friendIter.next());
			}
			uncoveredNodes.add(i);
		}

		Set<Integer> coveredNodes = new HashSet<Integer>();
		int i = 0;

		// while there is still an uncovered node in the graph
			while (!uncoveredNodes.isEmpty()){
				// find the nodes with the highest number of edges
					int v = findMaxNode(uncoveredNodes);
					// add node to dominating set
					dominatingSet.add(v);
					//Mark all the nodes v is connected to as covered
					removeEdge(v, uncoveredNodes, coveredNodes);
					//display progress
					i++;
					if (i > 1000){
						System.out.println(" Nodes left => " + uncoveredNodes.size());
						i = 0;
					}

			}

		return dominatingSet;
	}

	/*
	 * Project week 3 - assignment easier question: To find the dominating set of nodes in a network
	 * find the node with the maximum number of edges
	 */
	public Integer findMaxNode(Set<Integer> uncoveredNodes){
		int maxNode = 0;
		int countEdges = 0;
		Iterator<Integer> uncovIter = uncoveredNodes.iterator();
		while (uncovIter.hasNext()){
			int i = uncovIter.next();
			Set friends = adjListsMap.get(i);
			if (friends!=null){
				if (friends.size() > countEdges){
					countEdges = friends.size();
					maxNode = i;
				}
			} else
				{
				//ignore this node
				uncovIter.remove();
				}
		}

		return maxNode;
	}

	/*
	 * Project week 3 - assignment easier question: To find the dominating set of nodes in a network
	 * remove all the edges connected to this node
	 */
	public void removeEdge(Integer v, Set<Integer> uncoveredNodes, Set<Integer> coveredNodes){
		Iterator<Integer> graphIter = uncoveredNodes.iterator();

		// Remove the node v from the edges and mark the node on other end of v as covered
		while (graphIter.hasNext()){
			int i = graphIter.next();
			HashSet<Integer> friends = adjListsMap.get(i);
			if (friends!= null){
				Iterator<Integer> friendIter = friends.iterator();
				while (friendIter.hasNext()){
					int f = friendIter.next();
					if (f == v){
						friendIter.remove();
						graphIter.remove();
						if (!coveredNodes.contains(i)){
							coveredNodes.add(i);
						}
					}

				}
			} else
					{
						graphIter.remove();
						if (!coveredNodes.contains(i)){
							coveredNodes.add(i);
						}
					}
		}

		//Mark all the edges connected to v as covered and remove it from the uncovered nodes
		HashSet<Integer> friendsOfv = adjListsMap.get(v);
		if (friendsOfv != null){
			Iterator<Integer> friendIter = friendsOfv.iterator();
			while (friendIter.hasNext()){
				int f = friendIter.next();
				friendIter.remove();
				if (uncoveredNodes.contains(f)){
					uncoveredNodes.remove(f);
				}
				if (!coveredNodes.contains(f)){
					coveredNodes.add(f);
				}
			}
		}
		//Remove v from uncovered nodes
		if (uncoveredNodes.contains(v)){
			uncoveredNodes.remove(v);
		}
		// Add node v to uncovered nodes
		adjListsMap.remove(v);
		//Mark v as covered
		if (!coveredNodes.contains(v)){
			coveredNodes.add(v);
		}


	}

	/*
	 * Capstone project: Harder question - To detect sub-communities in a social network
	 * @param : Integer value of the highest betweenness score of an edge in the graph
	 * Method will iterate till the acceptable betweenness score is achieved
	 * For each iteration it will remove the edges that have the highest betweenness score
	 */

	public void detectCommunity(int maxAllowed){
		int i = 0;
		//perform the community detection as per the max betweenness score allowed
		while (maxBetweenness > (long)maxAllowed){
			maxBetweenness = (long) 00000;
			HashSet<Edge> maxBetweennessEdge = new HashSet<Edge>();
			//Get the betweenness score of all the edges in the graph
			HashMap<Edge, Double> betweenness = getBetweenness();
			//find the max betweenness score in the graph
			for (Edge e: betweenness.keySet()){
				long betweenScore = Math.round(betweenness.get(e));
				if ( betweenScore >= maxBetweenness){
					maxBetweenness = betweenScore;
				}
			}
			//Check if the max betweenness is below the allowed limit, if so exit the loop
			if (maxBetweenness <= (long)maxAllowed){
				break;
			}
			//debug
			System.out.println("\n Iteration = " + i++ + "\t Highest betweenness score = " + maxBetweenness);
				//find the edges with the highest betweenness score and add it to a Set
				for (Edge e1: betweenness.keySet()){
					if (Math.round(betweenness.get(e1)) == maxBetweenness){
						maxBetweennessEdge.add(e1);
					}
				}

				//remove the edges with the highest betweenness score from the graph
				Iterator<Edge> maxBetweenIter = maxBetweennessEdge.iterator();
				while (maxBetweenIter.hasNext()){
					Edge e2 = maxBetweenIter.next();
					edges.remove(e2);
					e2.getStartNode().getEdges().remove(e2);
				}

			} //end while loop
	}


	/*
	 * Capstone project : Harder question: Calculate the betweenness of all the edges of the graph
	 * We will use the brandes algorithm for finding the betweenness centrality of each node modified for edges
	 */
	public HashMap<Edge, Double> getBetweenness(){
		HashMap<Edge, Double> betweenScore = new HashMap<Edge, Double>();
		//Initialize the betweenness score to zero for all the edges
		for (Edge e : getEdgeList()){
			betweenScore.put(e,0.0);
		}

		for (Node startNode : getNodeList()){

			// define Stack for tracing the nodes
			Stack<Node> traceStack = new Stack<Node>();
			//define visited set
			HashSet<Node> visited = new HashSet<Node>();

			//store the edges in the shortest path to a node
			HashMap<Node,LinkedList<Edge>> shortestPEdges = new HashMap<Node,LinkedList<Edge>>();

			// store the number of shortest paths to a node
			HashMap<Node,Integer> numShortestPaths = new HashMap<Node,Integer>();

			// store length of shortest path to a node
			HashMap<Node,Double> shortestLength = new HashMap<Node,Double>();

			//Initialize all the HashMaps

			for (Node iNode: getNodeList()){
				numShortestPaths.put(iNode,0);
				shortestLength.put(iNode,Double.POSITIVE_INFINITY);
				shortestPEdges.put(iNode,new LinkedList<Edge>());
			}

			//Start node will have 1 shortest path and 0 shortest length
			numShortestPaths.put(startNode,1);
			shortestLength.put(startNode, 0.0);
			visited.add(startNode);
			// Step 1 : Find the shortest path count from startNode
			// We will now search from start to all other nodes using bfs
			//create an empty queue
			Queue<Node> bfs = new LinkedList<Node>();
			bfs.add(startNode);

			while (!bfs.isEmpty()){
				Node current = bfs.remove();
				traceStack.push(current);

				// for all the edges of the node
				for (Edge e : current.getEdges()){

						Node friend = e.getEndNode();
						if (!visited.contains(friend)){
							visited.add(friend);
							bfs.add(friend);
						}
							//if friend is found for the first time add 1 to the path length of predecessor
								if (shortestLength.get(friend) > shortestLength.get(current) + 1){
									shortestLength.put(friend,shortestLength.get(current)+1);
									numShortestPaths.put(friend,numShortestPaths.get(current));
									shortestPEdges.get(friend).add(e);
								} else
									{
									//if shortest path found to friend
										if (shortestLength.get(friend) == shortestLength.get(current)+1){
											numShortestPaths.put(friend,numShortestPaths.get(current)+
																		numShortestPaths.get(friend));
											shortestPEdges.get(friend).add(e);
											}
									}
				}
			}   // end bfs while loop

			// Step 2: Accumulate the dependencies visiting the nodes passed thru (Stack in bfs)
			//initialize the dependency HashMap for all edges
			HashMap<Edge, Double> dep = new HashMap<Edge, Double>();
			for (Edge e: getEdgeList()){
				dep.put(e,0.0);
			}


			while (!traceStack.isEmpty()){
				Node w = traceStack.pop();
				Double depNum = 0.0;
				for(Edge e : w.getEdges()){
						depNum += dep.get(e);
				}

				Iterator<Edge> shortestPIter = shortestPEdges.get(w).listIterator(0);
				while (shortestPIter.hasNext()){
						Edge e = shortestPIter.next();
						dep.put(e,((double)numShortestPaths.get(e.getStartNode())
								/(double)numShortestPaths.get(w))*(1+depNum));
						if (!w.equals(startNode)){
							betweenScore.put(e,betweenScore.get(e)+dep.get(e));
						}
				}
			} // end accumulation while loop
		} //end outer for loop

		return betweenScore;
	}
	/*
	 * main method - used for testing
	 */
	public static void main(String[] args)
	{
		System.out.print("Making a new Graph...");
		CapGraph theGraph = new CapGraph();
		System.out.print("DONE. \nLoading the Graph...");
		GraphLoader.loadGraph(theGraph,"data/small_test8_graph.txt");
		System.out.println(" (size " + theGraph.getNumVertices() + "+" + theGraph.getNumEdges() + " Graph nodes and edges )");
//		System.out.println(theGraph.printGraph());
		System.out.println("\n Printing graph: ");
		for (Node i: theGraph.getNodeList()){
			System.out.println("\n\t" + i.toString());
			for (Edge e: i.getEdges()){
				System.out.print("\t" + e.toString());
			}
		}
		System.out.println("\n End printing graph");


/*		Graph TGraph = new CapGraph();
		TGraph = theGraph.transposeGraph();
		System.out.println(TGraph.printGraph());
		System.out.println("End printing Transpose graph");
*/
/*		Set<Integer> dominatingSet = theGraph.getDominatingSet();
		System.out.println("dominating set size = " + dominatingSet.size());
		Iterator<Integer> dIter = dominatingSet.iterator();
		int i = 0;
		while (dIter.hasNext()){
			System.out.print("\t" + dIter.next());
			i++;
			if (i > 20){
				System.out.println("\n");
				i = 0;
			}
		}
*/
/*		HashMap<Edge, Double> betweenNess = theGraph.getBetweenness();
		for (Edge e : betweenNess.keySet()){
			System.out.println(e.toString() + " betweenness: " + betweenNess.get(e));
		}
*/
		theGraph.detectCommunity(10);

		System.out.println("Printing graph again: ");
		System.out.println(" (size " + theGraph.getNumVertices() + "+" + theGraph.getNumEdges() + " Graph nodes and edges )");

		for (Node i: theGraph.getNodeList()){
			System.out.println("\n\t" + i.toString());
			for (Edge e: i.getEdges()){
				System.out.print("\t" + e.toString());
			}
		}
		System.out.println("\n End printing graph again");
		System.out.println("\n" + "DONE.");


	} //end main


}
