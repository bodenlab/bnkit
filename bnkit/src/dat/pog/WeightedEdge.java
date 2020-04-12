package dat.pog;

public abstract class WeightedEdge implements Edge {
    private double weight = 0;

    public void setWeight(double weight) {
        this.weight = weight;
    }

    public double getWeight() {
        return weight;
    }
}
