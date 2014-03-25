/*
 * Copyright (C) 2014 M. Boden et al.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bn;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;

/**
 * Mixture of probability distributions.
 * 
 * @author mikael
 */
public class MixtureDistrib implements Distrib {

    final Map<Distrib, Double> mixture;
    private double density;
    Random rand = new Random();
    
    public MixtureDistrib(Distrib d1, double weight1) {
        mixture = new HashMap<>();
        mixture.put(d1, weight1);
        density = weight1;
    }
    
    public MixtureDistrib(Distrib d1, double weight1, Distrib d2, double weight2) {
        mixture = new HashMap<>();
        mixture.put(d1, weight1);
        mixture.put(d2, weight2);
        density = weight1 + weight2;
    }
    
    public double addDistrib(Distrib d2, double weight2) {
        Double prev_weight = mixture.get(d2);
        if (prev_weight == null)
            mixture.put(d2, weight2);
        else
            mixture.put(d2, prev_weight.doubleValue() + weight2);
        density += weight2;
        return density;
    }
    
    @Override
    public double get(Object value) {
        double p = 0.0;
        for (Map.Entry<Distrib, Double> entry : mixture.entrySet()) {
            p += entry.getKey().get(value) * entry.getValue().doubleValue();
        }
        return p;
    }

    @Override
    public Object sample() {
        double y = rand.nextDouble() * density;
        Distrib current = null;
        double p = 0.0;
        for (Map.Entry<Distrib, Double> entry : mixture.entrySet()) {
            current = entry.getKey();
            p += entry.getValue().doubleValue();
            if (p >= y)
                break;
        }
        if (current == null)
            throw new RuntimeException("Invalid MixtureDistrib");
        return current.sample();
    }
    
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("^" + mixture.size());
        for (Map.Entry<Distrib, Double> entry : mixture.entrySet())
            sb.append("{" + entry.getKey() + "*" + String.format("%3.1e", entry.getValue()) + "}");
        return sb.toString();
    }
    
}
