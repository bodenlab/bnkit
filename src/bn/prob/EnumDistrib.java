/*
    bnkit -- software for building and using Bayesian networks
    Copyright (C) 2014  M. Boden et al.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package bn.prob;

import bn.Distrib;
import dat.Domain;
import dat.Enumerable;
import json.JSONArray;
import json.JSONException;
import json.JSONObject;

import java.util.Map;
import java.util.Random;

/**
 * A probability distribution for enumerable variables.
 *
 * @author mikael
 */
public class EnumDistrib implements Distrib, Domain {

    private final double[] distrib;
    private final Enumerable domain;
    private boolean normalised = false;
    private boolean valid = false;
    private Random rand = new Random();

    /**
     * Initialise a probability distribution for a specified (enumerable)
     * domain. Does not assign any probabilities (only 0) so cannot be used
     * until at least one value of the domain has been associated with a
     * probability.
     *
     * @param discrete the categorical domain (that defines enumerable values)
     */
    public EnumDistrib(Enumerable discrete) {
        domain = discrete;
        distrib = new double[discrete.size()];
    }

    /**
     * Initialise a probability distribution for a specified (enumerable)
     * domain.
     *
     * @param discrete the categorical domain (that defines enumerable values)
     * @param distrib the probability distribution, or values proportional to
     * the same (will normalise)
     */
    public EnumDistrib(Enumerable discrete, double... distrib) {
        domain = discrete;
        if (distrib == null) {
            this.distrib = new double[discrete.size()];
        } else if (distrib.length == 0) {
            this.distrib = new double[discrete.size()];
        } else {
            this.distrib = new double[distrib.length];
            for (int i = 0; i < distrib.length; i++) {
                this.distrib[i] = distrib[i];
            }
            valid = true;
            normalise();
        }
    }

    /**
     * Initialise a probability distribution for a specified (enumerable)
     * domain.
     *
     * @param discrete the categorical domain (that defines enumerable values)
     * @param count counts of values (will normalise)
     */
    public EnumDistrib(Enumerable discrete, int... count) {
        domain = discrete;
        if (count.length == discrete.size()) {
            this.distrib = new double[discrete.size()];
            for (int i = 0; i < distrib.length; i++) {
                this.distrib[i] = count[i];
            }
            valid = true;
            normalise();
        } else
            throw new RuntimeException("Invalid EnumDistrib");
    }

    /**
     * Initialise a probability distribution from a map that both defines the
     * values of a domain (which will be constructed internally) and the
     * corresponding probabilities.
     *
     * @param enumdistrib map with values and their probabilities (will be
     * normalised if necessary)
     */
    public EnumDistrib(Map<Object, Double> enumdistrib) {
        Object[] mydomain = new Object[enumdistrib.size()];
        double[] mydistrib = new double[enumdistrib.size()];
        int i = 0;
        for (Map.Entry<Object, Double> entry : enumdistrib.entrySet()) {
            mydomain[i] = entry.getKey();
            mydistrib[i] = entry.getValue().doubleValue();
            i++;
        }
        this.domain = new Enumerable(mydomain);
        this.distrib = mydistrib;
        valid = true;
        normalise();
    }
    
    /**
     * Initialise a probability distribution from a map that both defines the
     * values of a domain (which will be constructed internally) and the
     * corresponding probabilities.
     *
     * @param enumdistrib map with values and their probabilities (will be
     * normalised if necessary)
     * @param discrete the categorical domain (that defines enumerable values)
     */
    public EnumDistrib(Map<Object, Double> enumdistrib, Enumerable discrete) {
        double[] mydistrib = new double[enumdistrib.size()];
        int i = 0;
        for (Map.Entry<Object, Double> entry : enumdistrib.entrySet()) {
            mydistrib[i] = entry.getValue().doubleValue();
            i++;
        }
        this.domain = discrete;
        this.distrib = mydistrib;
        valid = true;
        normalise();
    }
    
    /**
     * Get domain of a conjugate distribution.
     * */
    public Enumerable getDomain() {
        return domain;
    }
    /**
     * Set the seed for the random number generation that may be used for
     * sampling etc.
     *
     * @param seed
     */
    public void setSeed(long seed) {
        rand = new Random(seed);
    }

    /**
     * Create a JSON representation of distribution
     * @return
     */
    public JSONObject toJSON() {
        JSONObject json = new JSONObject();
        json.put("Domain", domain.toJSON());
        json.put("Pr", new JSONArray(this.get()));
        return json;
    }

    public static EnumDistrib fromJSON(JSONObject json) {
        Enumerable domain = Enumerable.fromJSON(json.getJSONObject("Domain"));
        double[] pr = new double[domain.size()];
        JSONArray prarr = json.getJSONArray("Pr");
        assert(pr.length == prarr.length());
        for (int i = 0; i < prarr.length(); i ++)
            pr[i] = prarr.getDouble(i);
        return new EnumDistrib(domain, pr);
    }

    public String toString() {
        if (!isNormalised()) {
            normalise();
        }
        StringBuffer sbuf = new StringBuffer("<");
        for (int i = 0; i < distrib.length; i++) {
            sbuf.append(String.format("%4.2f ", distrib[i]));
        }
        sbuf.replace(sbuf.length() - 1, sbuf.length() - 1, ">");
        return sbuf.toString();
    }

    public static EnumDistrib parseEnumDistrib(String str, Enumerable dom) {
        String[] parts = str.split(",");
        if (parts.length != dom.size())
            throw new RuntimeException("Invalid specification of distribution: " + str);
        double[] values = new double[parts.length];
        
        for (int i = 0; i < values.length; i ++) {
            try {
                values[i] = Double.valueOf(parts[i]);
            } catch (NumberFormatException e) {
                System.err.println("Invalid parameter for distribution: " + parts[i] + " in " + str);
            }
        }
        return new EnumDistrib(dom, values);
    }
    
    /**
     * Retrieve the probability of the specified value
     *
     * @param value the instance for which a probability is sought
     * @return the probability value
     */
    public double get(Object value) {
        if (!isNormalised()) {
            normalise();
        }
        int index = domain.getIndex(value);
        return distrib[index];
    }

    /**
     * Retrieve the probability associated with the value at the specified index
     * of the domain.
     *
     * @param index of value
     * @return the probability value
     */
    public double get(int index) {
        if (!isNormalised()) {
            normalise();
        }
        return distrib[index];
    }

    public double[] get() {
        if (!isNormalised()) {
            normalise();
        }
        return distrib;
    }

    /**
     * Get the index of the value that is assigned the highest probability
     * @return the index
     */
    public int getMaxIndex() {
        int max_index = 0;
        for (int i = 1; i < distrib.length; i++) {
            if (distrib[max_index] < distrib[i])
                max_index = i;
        }
        return max_index;
    }

    /**
     * Get the value that is assigned the highest probability
     * @return the value
     */
    public Object getMax() {
        int max_index = getMaxIndex();
        return domain.get(max_index);
    }

    /**
     * Choose a value from the domain, with a probability based on the
     * distribution.
     *
     * @return a value of the domain
     */
    public Object sample() {
        if (!isNormalised()) {
            normalise();
        }
        double p = rand.nextDouble();
        double sum = 0;
        for (int i = 0; i < distrib.length; i++) {
            sum += distrib[i];
            if (sum >= p) {
                return this.domain.get(i);
            }
        }
        throw new RuntimeException("Invalid sampling");
    }

    /**
     * Set the probability of one specific value in the domain. The distribution
     * is flagged as not normalised as a result. If a series of values are set
     * which do not add to 1, the distribution needs to be normalised manually.
     *
     * @param value the value of the domain which is assigned a probability
     * @param prob the probability
     */
    public void set(Object value, double prob) {
    	if (!isValid()) {
    		for (int i = 0; i < distrib.length; i ++)
    			distrib[i] = 0;
    		setValid(true);
    	}
        normalised = false;
        int index = domain.getIndex(value);
        distrib[index] = prob;
    }

    /**
     * Set the probabilities associated with this enumerable domain. If the
     * probabilities do not add up to 1, they will be normalised.
     *
     * @param distrib the probability distribution, or values proportional to
     * their probabilities.
     */
    public void set(double[] distrib) {
        if (distrib.length == this.distrib.length) {
            for (int i = 0; i < distrib.length; i++) {
                this.distrib[i] = distrib[i];
            }
            normalise();
            setValid(true);
        } else {
            throw new RuntimeException("Invalid size distribution");
        }
    }

    public boolean isNormalised() {
        return normalised;
    }

    public boolean isValid() {
        return valid;
    }
    
    public void setValid(boolean status) {
    	valid = status;
    }

    public void normalise() {
        double sum = 0.0;
        for (int i = 0; i < distrib.length; i++) {
            sum += distrib[i];
        }
        for (int i = 0; i < distrib.length; i++) {
            distrib[i] /= sum;
        }
        normalised = true;
    }

    public static EnumDistrib uniform(Enumerable discrete) {
        EnumDistrib d = new EnumDistrib(discrete);
        for (int i = 0; i < d.distrib.length; i++) {
            d.distrib[i] = 1.0 / discrete.size();
        }
        d.normalised = true;
        return d;
    }

    public static EnumDistrib random(Enumerable discrete, long seed) {
        Random rand = new Random(seed);
        EnumDistrib d = new EnumDistrib(discrete);
        double sum = 0.0;
        for (int i = 0; i < d.distrib.length; i++) {
            d.distrib[i] = rand.nextDouble();
            sum += d.distrib[i];
        }
        for (int i = 0; i < d.distrib.length; i++) {
            d.distrib[i] /= sum;
        }
        d.normalised = true;
        return d;
    }

    public static EnumDistrib random(Enumerable discrete) {
        Random rand = new Random();
        return random(discrete, rand.nextInt());
    }

    public double getEntropy() {
        if (!isNormalised()) 
            normalise();
        double sum = 0;
        for (double p : distrib)
            sum += p * Math.log(p);
        return -sum;
    }
    
    @Override
    public boolean isValid(Object value) {
        try {
            EnumDistrib x = (EnumDistrib) value;
            return true;
        } catch (ClassCastException e) {
            return false;
        }
    }

    public static double[] log2Prob(double[] logprobs) {
        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        double[] p = new double[logprobs.length];
        for (int i = 0; i < logprobs.length; i ++) {
            if (logprobs[i] < min) 
                min = logprobs[i];
            if (logprobs[i] > max)
                max = logprobs[i];
        }
        double sum = 0;
        for (int i = 0; i < logprobs.length; i ++) {
            p[i] = Math.exp(logprobs[i] - max);
            sum += p[i];
        }        
        for (int i = 0; i < logprobs.length; i ++)
            p[i] /= sum;
        return p;
    }
}
