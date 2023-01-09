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

import json.JSONArray;
import json.JSONObject;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * This class represents a (countable) domain of values that discrete variables
 * can take.
 *
 * @author mikael
 */
public class Enumerable implements Domain {

    final int order;
    final Object[] values;
    final Class datatype;

    public Enumerable(int order) {
        if (order < 2) {
            throw new RuntimeException("An Enumerable must have at least two values");
        }
        this.order = order;
        this.values = new Integer[order];
        for (int i = 0; i < order; i++) {
            this.values[i] = i;
        }
        this.datatype = Integer.class;
    }

    /**
         * Check values to what class they all belong
     * @return the class to which all values belong
     * @throws RuntimeException if values contain more than one class
     * */

    public Enumerable(Object[] values) {
        this.values = values;
        this.order = values.length;
        if (order < 2) {
            throw new RuntimeException("An Enumerable must have at least two values");
        }
        Class c = null;
        for (Object v : values) {
            if (c == null)
                c = v.getClass();
            else if (c != v.getClass())
                throw new RuntimeException("Invalid mix of values " + c.toString() + " and " + v.getClass().toString());
        }
        this.datatype = c;
    }

    /**
     */
    private Class getDatatype() {
        return this.datatype;
    }

    public int size() {
        return order;
    }
    
    public boolean equals(Enumerable enumerable) {
    	if (size() != enumerable.size()) {
    		return false;
    	}
        for (Object sym : enumerable.getValues()) {
            if (!this.isValid(sym))
                return false;
    	}
    	return true;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        for (int i = 0; i < order; i++)
            hash = 31 * hash + this.values[i].hashCode();
        return hash;
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
    public static Enumerable gap_ext = new Enumerable(new Character[]{'G','C','?'});

    public static Map<String, Enumerable> predef = new HashMap<>();
    public static Map<Enumerable, String> predef_reverse = new HashMap<>();

    // Instantiating the static map
    static
    {
        predef.put("Boolean", bool);
        predef.put("DNA", nacid);
        predef.put("RNA", nacidRNA);
        predef.put("DNA with N", nacidwn);
        predef.put("RNA with N", nacidwnRNA);
        predef.put("Protein", aacid);
        predef.put("Protein with X", aacidwx);
        predef.put("Protein with gap", aacid_ext);
        predef_reverse.put(bool, "Boolean");
        predef_reverse.put(nacid, "DNA");
        predef_reverse.put(nacidRNA, "RNA");
        predef_reverse.put(nacidwn, "DNA with N");
        predef_reverse.put(nacidwnRNA, "RNA with N");
        predef_reverse.put(aacid, "Protein");
        predef_reverse.put(aacidwx, "Protein with X");
        predef_reverse.put(aacid_ext, "Protein with gap");
    }

    public JSONObject toJSON() {
        JSONObject jobj = new JSONObject();
        String defID = predef_reverse.get(this);
        if (defID != null) { // predefined
            jobj.put("Predef", defID);
            return jobj;
        }
        // not predefined, so provide details
        jobj.put("Size", order);
        jobj.put("Values", new JSONArray(values));
        jobj.put("Datatype", datatype.getSimpleName());
        return jobj;
    }

    public static Enumerable fromJSON(JSONObject jobj) {
        String defID = jobj.optString("Predef", null);
        if (defID != null) {
            Enumerable match = predef.get(defID);
            if (match == null) {
                throw new RuntimeException("Invalid predefined domain name: \"" + defID + "\". Please check dat.Enumerable");
            }
            return match;
        }
        JSONArray varr = jobj.getJSONArray("Values");
        //Class dtype = (Class)jobj.get("Datatype");
        Object[] vals = new Object[varr.length()];
        for (int i = 0; i < varr.length(); i++)
            vals[i] = varr.get(i);
        return new Enumerable(vals);
    }

    /**
     * A class to associate an enumerable with a value.
     * Useful for passing enumerable/value pairs to methods where multiples of these are required.
     */
    public static class Designation {
        public final Enumerable enumerable;
        public final Object val;
        public Designation(Enumerable enumerable, Object val) {
            this.enumerable = enumerable;
            this.val = val;
        }
        public static Designation[] array(Enumerable[] enumerables, Object[] vals) {
            Designation[] ret = new Designation[enumerables.length];
            for (int i = 0; i < enumerables.length && i < vals.length; i ++)
                ret[i] = new Designation(enumerables[i], vals[i]);
            return ret;
        }
        public static <T extends Enumerable> Designation[] array(List<T> enumerables, Object[] vals) {
            Designation[] ret = new Designation[enumerables.size()];
            for (int i = 0; i < ret.length && i < vals.length; i ++)
                ret[i] = new Designation(enumerables.get(i), vals[i]);
            return ret;
        }
        public static Map<Enumerable, Object> toMap(Designation[] assign_array) {
            Map<Enumerable, Object> map = new HashMap<>();
            for (Designation assign : assign_array)
                map.put(assign.enumerable, assign.val);
            return map;
        }
        public static Map<Enumerable, Object> toMap(Collection<Designation> assign_list) {
            Map<Enumerable, Object> map = new HashMap<>();
            for (Designation assign : assign_list)
                map.put(assign.enumerable, assign.val);
            return map;
        }
        public String toString() {
            return enumerable.toString() + "=" + val.toString();
        }
    }

    /**
     * A factory method for creating a single variable/value pair.
     * @param enumerable
     * @param val
     * @return the variable paired with a value
     */
    public static Designation assign(Enumerable enumerable, Object val) {
        return new Designation(enumerable, val);
    }


}
