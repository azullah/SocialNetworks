package graph;

import java.util.HashSet;
import java.util.Set;

public class Node {

	private int value;

	private HashSet<Edge> edges;

	public Node(int v){
		edges = new HashSet<Edge>();
		value = v;
	}

	public int getValue(){
		   return value;
	}

	void addEdge(Edge e)
	{
		edges.add(e);
	}

	void setValue(int s){
		value = s;
	}

	public Set<Node> getNeighbors()
	{
		Set<Node> neighbors = new HashSet<Node>();
		for (Edge edge : edges) {
			neighbors.add(edge.getOtherNode(this));
		}
		return neighbors;
	}

	public HashSet<Edge> getEdges(){
		return edges;
	}

	public String toString()
	{
		String toReturn = "Node: ";
		toReturn += "\t" + value;
		return toReturn;
	}

	public boolean equals(Object o)
	{
		if (!(o instanceof Node) || (o == null)) {
			return false;
		}
		Node node = (Node)o;
		return (node.getValue() == this.value);
	}

	public int hashCode()
	{
		return value;
	}

}
