/*
 * bnkit -- software for building and using Bayesian networks
 * Copyright (C) 2014  M. Boden et al.
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
package dat;

import json.JSONArray;
import json.JSONObject;

/**
 *
 * @author mikael
 * @param <E> the domain that each element belongs to
 */
public class SeqDomain<E extends Domain> implements Domain {
    
    protected final E elementType;
    protected Object[] arr = null;
    
    public SeqDomain(E elementType) {
        this.elementType = elementType;
    }

    @Override
    public boolean isValid(Object value) {
        try {
            Object[] arr = (Object[]) value;
            for (int i = 0; i < arr.length; i ++) {
                if (!elementType.isValid(arr[i]))
                    return false;
            }
        } catch (ClassCastException e) {
            return false;
        }
        return true;
    }
    
    public boolean set(Object[] value) {
        if (this.isValid(value)) {
            //System.err.println("Invalid character found: " + value);
            this.arr = value;
            return true;
        }
        return false;
    }
    
    public Object[] get() {
        return arr;
    }

    public Object get(int idx) {
        return arr[idx];
    }

    @Override
    public boolean equals(Object other) {
        try {
            SeqDomain odom = (SeqDomain) other;
            if (odom.arr.length != this.arr.length)
                return false;
            if (!odom.getType().equals(this.getType()))
                return false;
            for (int i = 0; i < arr.length; i ++)
                if (arr[i] != null) {
                    if (!arr[i].equals(odom.arr[i]))
                        return false;
                } else if (odom.arr[i] != null) // this null, other is NOT
                    return false;
        } catch (ClassCastException e) {
            return false;
        }
        return true;
    }

    public int length() {
        if (arr != null)
            return arr.length;
        else
            return 0;
    }
    
    public E getType() {
        return elementType;
    }
    
    public static SeqDomain<Continuous> cont_seq = new SeqDomain<>(new Continuous());

    public static void main(String[] args) {
        SeqDomain seq = new SeqDomain<>(Enumerable.aacid);
        System.out.println(seq.isValid(new Character[] {'A','A','A','G','G','Q','T','C','C','C','A','C','T'}));
        
    }

    /**
     * Count the number of null elements in the sequence.
     *
     * @return the number of null elements in the sequence
     */
    public int getNNulls() {
        int cnt = 0;
        for (int i = 0; i < arr.length; i ++) {
            if (arr[i] == null)
                cnt += 1;
        }
        return cnt;
    }

    /**
     * Counts the number of null sequences that start in the array.
     *
     * This method iterates through the array and counts the number of times
     * a null sequence starts. A null sequence is considered started when a
     * null element is encountered.
     *
     * @return the number of null sequences that start in the array
     */
    public int getNNullStarts() {
        int cnt = 0;
        boolean inGap = false;
        for (int i = 0; i < arr.length; i ++) {
            if (arr[i] == null && !inGap) {
                cnt += 1;
                inGap = true;
            } else
                inGap = false;
        }
        return cnt;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (Object val : get())
            sb.append("|"+val.toString());
        sb.append("|");
        return sb.toString();
    }
    @Override
    public JSONObject toJSON() {
        JSONObject jobj = new JSONObject();
        jobj.put("Predef", "Sequence");
        return jobj;
    }

}
