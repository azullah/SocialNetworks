package graph;

// Capstone project: New class created for storing edges

public class Edge {
	private Node start;
	private Node end;

	public Edge (Node st, Node en){
		start = st;
		end = en;
	}

	public Node getStartNode(){
		return start;
	}

	public Node getEndNode(){
		return end;
	}

	public Node getOtherNode(Node node)
	{
		if (node.equals(start))
			return end;
		else if (node.equals(end))
			return start;
		throw new IllegalArgumentException("Looking for " +
			"a node that is not in the edge");
	}
	public String toString()
	{
		String toReturn = "EDGE between ";
		toReturn += "\t" + start.getValue();
		toReturn += "\t" + end.getValue();
		return toReturn;
	}

}
