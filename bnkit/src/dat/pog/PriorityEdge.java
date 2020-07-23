package dat.pog;

public class PriorityEdge implements WeightedEdge {

    private double weight = 1;  // default
    public PriorityEdge() {

    }

    @Override
    public double getWeight() {
        return weight;
    }

}
