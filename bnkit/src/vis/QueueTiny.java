/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vis;

/**
 *
 * @author ariane
 */
public class QueueTiny<Data> {

    /**
     * The total cost so far to reach the state.
     */
    public double totalCost;
    /**
     * Current position of the person on the map
     */
    private Integer currentPosition;
    /**
     * Current position of the person on the map
     */
    private Integer previousPosition;
    /**
     * A heuristic estimate of the remaining cost to the goal; this will be zero
     * if no heuristic is being used
     */
    public double heuristicEstimate;
    /**
     * The id. Allows us to know if it is the same state as before
     */
    public int id;

    /**
     * Constructs a queue entry with the given parameters.
     *
     * @param id
     * @param totalCost the total cost so far to reach the state.
     * @param heuristicEstimate an estimate of the cost to the goal; this should
     * be zero if no heuristic is used.
     * @param currentPosition
     * @param previousPosition
     */
    public QueueTiny(int id, double totalCost, double heuristicEstimate, Integer currentPosition, Integer previousPosition) {
        this.totalCost = totalCost;
        this.heuristicEstimate = heuristicEstimate;
        this.currentPosition = currentPosition;
        this.previousPosition = previousPosition;
        this.id = id;
    }

    /**
     * Returns the id.
     *
     * @return the id of the queue entry.
     */
    public int getId() {
        return id;
    }

    /**
     * Returns the total cost so far.
     *
     * @return the total cost so far.
     */
    public double getTotalCost() {
        return totalCost;
    }

    /**
     * Returns a heuristic estimate of the remaining cost (zero if N/A).
     *
     * @return a heuristic estimate of the remaining cost (zero if N/A).
     */
    public double getHeuristicEstimate() {
        return heuristicEstimate;
    }

    /**
     * Returns the current position.
     *
     * @return
     */
    public Integer getPreviousPosition() {
        return previousPosition;
    }

    /**
     * Returns the current position.
     *
     * @return
     */
    public Integer getCurrentPosition() {
        return currentPosition;
    }

    /**
     * Implements a comparison between queue entries, based on the sum of the
     * total cost so far and the heuristic estimate. This makes A* and UCS
     * function naturally with a PriorityQueue.
     * @param arg0
     */
    public int compareTo(QueueTiny<Data> arg0) {
        return Double.compare(this.totalCost + this.heuristicEstimate, arg0.totalCost + arg0.heuristicEstimate);
    }

}
