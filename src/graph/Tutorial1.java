package graph;
import org.graphstream.graph.*;
import org.graphstream.graph.implementations.*;

import util.GraphLoader;
import java.util.HashMap;
import org.graphstream.algorithm.ConnectedComponents;
import java.util.*;

class Tutorial1 {

   public static String styleSheet =
		   "node {size: 10px;fill-mode: dyn-plain;fill-color: black, red;z-index: 0;}" +
		   "edge {shape: line;fill-color: blue;arrow-size: 3px, 2px;}";
/*	            "node {" +
//	        	" 	text-mode: display;" +
	            "	fill-mode: dyn-plain;" +
	        	"	fill-color: white, black;" +
	        	"	size-mode: dyn-size;" +
	            "}" +
	            "node.marked {" +
	            "	fill-color: red;" +
	            "}";
*/

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		MultiGraph displayGraph = new MultiGraph("Tutorial 1");
		CapGraph theGraph = new CapGraph();

		System.out.print("\nLoading the Graph...");
		GraphLoader.loadGraph(theGraph,"data/facebook_2000.txt");
		System.out.print("\nDONE.");
		System.out.println(" (size " + theGraph.getNumVertices() + "+" + theGraph.getNumEdges() + " Graph nodes and edges )");

		System.out.print("\nDetecting communities...");
		theGraph.detectCommunity(500);
		HashMap<Edge, Double> betweenNess = theGraph.getBetweenness();

//		Set<Integer> dominatingSet = theGraph.getDominatingSet();
//		System.out.println("dominating set size = " + dominatingSet.size());

		for (Node i: theGraph.getNodeList()){
			String n = String.valueOf(i.getValue());
			displayGraph.addNode(n);

		}
		for (Edge e: theGraph.getEdgeList()){
			String s = String.valueOf(e.getStartNode().getValue()) + "-" + String.valueOf(e.getEndNode().getValue());
			displayGraph.addEdge(s, String.valueOf(e.getStartNode().getValue()),String.valueOf(e.getEndNode().getValue()));
//			displayGraph.getEdge(s).addAttribute("ui.label", String.valueOf(betweenNess.get(e)));
			displayGraph.getEdge(s).addAttribute("ui.color", 0);
		}

		   for (org.graphstream.graph.Node node : displayGraph.getEachNode()) {
//		        node.addAttribute("ui.label", " " + node.getId());
/*					if (dominatingSet.contains(Integer.parseInt(node.getId()))){
				        node.addAttribute("ui.color", 1);
//						System.out.println("dominating set node " + node.getId());
					} else */
				        node.addAttribute("ui.color", 1);
		    }

		   ConnectedComponents cc = new ConnectedComponents();
		   cc.init(displayGraph);

			System.out.print("\nNum of communities detected = " + cc.getConnectedComponentsCount());

		   displayGraph.addAttribute("ui.stylesheet", styleSheet);
		   displayGraph.addAttribute("ui.quality");
		   displayGraph.addAttribute("ui.antialias");
		   displayGraph.display();


	}

}
