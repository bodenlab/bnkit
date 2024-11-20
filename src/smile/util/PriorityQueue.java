/*
 * Copyright (c) 2010-2021 Haifeng Li. All rights reserved.
 *
 * Smile is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Smile is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Smile.  If not, see <https://www.gnu.org/licenses/>.
 */

package smile.util;

/**
 * Priority Queue for index items.
 *
 * @author Haifeng Li
 */
public class PriorityQueue {

    /**
     * The number of items in the queue.
     */
    private int n;
    /**
     * The d-ary heap or d-heap is a generalization of the binary heap data
     * structure whose non-leaf nodes have d children, instead of 2. Thus,
     * a binary heap is a 2-heap.
     */
    private final int d;
    /**
     * External array of priority.
     */
    private final double[] a;
    /**
     * The array of item indices.
     */
    private final int[] pq;
    /**
     * The inverse array qp allows the priority-queue to treat the array indices
     * as handles.
     */
    private final int[] qp;

    /**
     * Priority comparison of item i and j.
     * @param i item index
     * @param j item index
     */
    private boolean less(int i, int j) {
        return a[pq[i]] < a[pq[j]];
    }

    /**
     * Swap i and j items of pq and qp.
     * @param i item index
     * @param j item index
     */
    private void swap(int i, int j) {
        int t = pq[i];
        pq[i] = pq[j];
        pq[j] = t;
        qp[pq[i]] = i;
        qp[pq[j]] = j;
    }

    /**
     * fix up.
     */
    private void swim(int k) {
        while (k > 1 && less(k, (k + d - 2) / d)) {
            swap(k, (k + d - 2) / d);
            k = (k + d - 2) / d;
        }
    }

    /**
     * fix down.
     */
    private void sink(int k, int N) {
        int j;
        while ((j = d * (k - 1) + 2) <= N) {
            for (int i = j + 1; i < j + d && i <= N; i++) {
                if (less(i, j)) {
                    j = i;
                }
            }
            if (!(less(j, k))) {
                break;
            }
            swap(k, j);
            k = j;
        }
    }

    /**
     * Constructor. Default use a 3-heap.
     * @param a external array of priority. Lower value means higher priority.
     */
    public PriorityQueue(double[] a) {
        this(3, a);
    }

    /**
     * Constructor.
     * @param d d-heap.
     * @param a external array of priority. Lower value means higher priority.
     */
    public PriorityQueue(int d, double[] a) {
        this.d = d;
        this.a = a;
        this.n = 0;
        pq = new int[a.length + 1];
        qp = new int[a.length + 1];
    }

    /**
     * Returns true if the queue is empty.
     * @return true if the queue is empty.
     */
    public boolean isEmpty() {
        return n == 0;
    }

    /**
     * Insert a new item into queue.
     * @param v the index of item.
     */
    public void insert(int v) {
        pq[++n] = v;
        qp[v] = n;
        swim(n);
    }

    /**
     * Removes and returns the index of item with minimum value (highest priority).
     * @return the index of item with minimum value.
     */
    public int poll() {
        swap(1, n);
        sink(1, n - 1);
        return pq[n--];
    }

    /**
     * The value of item k is lower (higher priority) now.
     * @param k the item index.
     */
    public void lower(int k) {
        swim(qp[k]);
    }

    /**
     * The priority of item k has changed.
     * @param k the item index.
     */
    public void change(int k) {
        swim(qp[k]);
        sink(qp[k], n);
    }
}
