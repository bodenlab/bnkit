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
package dat;

/**
 * This class represents a (countable) domain of values that discrete variables
 * can take.
 *
 * @author mikael
 */
public class Enumerable implements Domain {

    final int order;
    final Object[] values;

    public Enumerable(int order) {
        if (order < 2) {
            throw new RuntimeException("An Enumerable must have at least two values");
        }
        this.order = order;
        this.values = new Integer[order];
        for (int i = 0; i < order; i++) {
            this.values[i] = i;
        }
    }

    public Enumerable(Object[] values) {
        this.values = values;
        this.order = values.length;
        if (order < 2) {
            throw new RuntimeException("An Enumerable must have at least two values");
        }
    }

    public int size() {
        return order;
    }
    
    public boolean equals(Enumerable enumerable) {
    	if(size() != enumerable.size()) {
    		return false;
    	}
    	
    	for(int i = 0; i < size(); i++) {
    		if(!get(i).equals(enumerable.get(i))) {
    			return false;
    		}
    	}
    	
    	return true;
    }

    /**
     * Retrieve the index of the value in the domain.
     * TODO: Improve speed for domains with many values.
     * @param value
     * @return 
     */
    public int getIndex(Object value) {
        if (value instanceof java.lang.String) {
            for (int i = 0; i < values.length; i++) {
                if (((String) value).equals(values[i])) {
                    return i;
                }
            }
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < values.length; i++)
                sb.append(values[i].toString() + ((i < values.length - 1)?", ":""));
            throw new RuntimeException("Value \"" + value.toString() + "\" unknown to enumerable domain " + this.toString() + " with values: " + sb.toString());
        } else {
            for (int i = 0; i < values.length; i++) {
                if (value.equals(values[i])) {
                    return i;
                }
            }
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < values.length; i++)
                sb.append(values[i].toString() + ((i < values.length - 1)?", ":""));
            throw new RuntimeException("Value \"" + value.toString() + "\" unknown to enumerable domain " + this.toString() + " with values: " + sb.toString());
        }
    }

    /**
     * Retrieve the index for each of the symbols in the array
     * @param value_array array of symbols
     * @return array of indices
     */
    public int[] getIndices(Object[] value_array) {
        int[] idx = new int[value_array.length];
        for (int i = 0; i < value_array.length; i ++) 
            idx[i] = getIndex(value_array[i]);
        return idx;
    }
    
    public Object get(int index) {
        return values[index];
    }

    public boolean isValid(Object value) {
        try {
            getIndex(value);
            return true;
        } catch (RuntimeException e) {
            return false;
        }
    }

    public Object[] getValues(){
        return this.values;
    }

    /**
     * Determine a unique hash key that can be used for encoding a word of enumerables.
     * @param arr the symbols making up the word
     * @return the key
     */
    public int getKey4Word(Object[] arr) {
        int sum1 = 0;
        int multiplier1 = 1;
        for (int i = 0; i < arr.length; i ++) {
            int idx = getIndex(arr[i]);
            sum1 += multiplier1 * idx;
            multiplier1 *= this.order;
        }
        return sum1;
    }

    /**
     * Determine a unique hash key that can be used for encoding a word of enumerables, here specified as a substring of an array.
     * @param arr the symbols making up the word
     * @param start the start index
     * @param end the end index
     * @return the key
     */
    public int getKey4Word(Object[] arr, int start, int end) {
        int sum1 = 0;
        int multiplier1 = 1;
        for (int i = start; i < arr.length && i < end; i ++) {
            int idx = getIndex(arr[i]);
            sum1 += multiplier1 * idx;
            multiplier1 *= this.order;
        }
        return sum1;
    }

    /**
     * Create the word (array of Object) that the specified key encodes.
     * @param key
     * @param k length of word
     * @return the word
     */
    public Object[] getWord4Key(int key, int k) {
        int maxmultip = (int)Math.pow(order, k);
        int multiplier = maxmultip / order;
        int remainder = key;
        Object[] word = new Object[k];
        for (int i = k - 1; i >= 0; i --) {
            int idx = remainder / multiplier;
            if (idx < 0 || idx >= order)
                throw new RuntimeException("Index invalid: " + key);
            word[i] = get(idx);
            remainder = remainder % multiplier;
            multiplier /= order;
        }
        return word;
    }

    public static String toString(Object[] word) {
        StringBuilder sb = new StringBuilder();
        for (Object sym : word)
            sb.append(sym.toString());
        return sb.toString();
    }
    public String toString(int key, int k) {
        Object[] word = getWord4Key(key, k);
        return Enumerable.toString(word);
    }
    
    public static Enumerable bool = new Enumerable(new Boolean[]{true, false});
    public static Enumerable nacid = new Enumerable(new Character[]{'A', 'C', 'G', 'T'});
    public static Enumerable nacidRNA = new Enumerable(new Character[]{'A', 'C', 'G', 'U'});
    public static Enumerable nacidwn = new Enumerable(new Character[]{'A', 'C', 'G', 'T', 'N'});
    public static Enumerable nacidwnRNA = new Enumerable(new Character[]{'A', 'C', 'G', 'U', 'N'});
    public static Enumerable secondaryStructure = new Enumerable(new Character[]{'H', 'C', 'E'});
    public static Enumerable aacid = new Enumerable(new Character[]{'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'});
    public static Enumerable aacidwx = new Enumerable(new Character[]{'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X'});
    public static Enumerable aacid_ext = new Enumerable(new Character[]{'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-'});
    public static Enumerable aacid_alt = new Enumerable(new Character[]{'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'});
    public static Enumerable gap_character = new Enumerable(new Character[]{'G','C'});

}
