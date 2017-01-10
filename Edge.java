//SURBHI LOHIA
//SL3893
//DATA STRUCTURES HOMEWORK 6
//MAY 3, 2016

public class Edge {

  public double distance;
  public Vertex source;
  public Vertex target;

  public Edge(Vertex vertex1, Vertex vertex2, double weight) {
    source = vertex1;
    target = vertex2; //adjacent vertex
    this.distance = weight;
  }

  public String toString() {
    return source + " - " + target;
  }
}