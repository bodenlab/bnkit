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

import dat.file.AlnReader;
import dat.file.FastaReader;
import json.JSONArray;
import json.JSONException;
import json.JSONObject;
import json.JSONString;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.Buffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 *
 * @author mikael
 * @param <E>
 */
public class EnumSeq<E extends dat.Enumerable> extends dat.SeqDomain<E> {

    String name;
    String info;
    
    public EnumSeq(E elementType) {
        super(elementType);
    }
    
    public void setName(String name) {
        this.name = name;
    }
    
    public void setInfo(String info) {
    	this.info = info;
    }
    
    public String getName() {
        return name;
    }

    public Object getFromIndex(int index) { return arr[index];}

    @Override
    public boolean equals(Object other) {
        try {
            EnumSeq oseq = (EnumSeq) other;
            if (name != null)
                if (!getName().equals(oseq.getName()))
                    return false;
            if (info != null)
                if (!getInfo().equals(oseq.getInfo()))
                    return false;
            return super.equals(other);
        } catch (ClassCastException e) {
            return false;
        }
    }

    @Override
    public E getType() {
        return super.getType();
    }
    
    public String getInfo() {
    	return info;
    }

    public JSONObject toJSONslim() {
        JSONObject json = new JSONObject();
        json.put("Name", getName());
        json.put("Seq", fromObjectArray(get()));
        return json;
    }
    public JSONObject toJSON() {
        JSONObject json = new JSONObject();
        json.put("Name", getName());
        json.put("Seq", fromObjectArray(get()));
        json.put("Datatype", getType().toJSON());
        return json;
    }

    public String toString() {
        StringBuffer sb = new StringBuffer();
        Object[] syms = get();
        for (int i = 0; i < syms.length; i ++) {
            if (syms[i] == null)
                sb.append('-');
            else
                sb.append(syms[i]);
        }
        return sb.toString();
    }

    public String toParsableString() {
        StringBuffer sb = new StringBuffer(length() > 0 ? "|" : "");
        Object[] syms = get();
        for (int i = 0; i < syms.length; i ++) {
            if (syms[i] == null)
                sb.append("null|");
            else
                sb.append(syms[i].toString() + "|");
        }
        return sb.toString();
    }

    public static JSONArray fromObjectArray(Object[] array) {
        JSONArray jsonarr = new JSONArray(array);
        return jsonarr;
    }

    /**
     * Converts a JSONArray of what-looks-like single-char strings to a character array, with null indicating the absence of string (i.e. gaps)
     * @param jsonArray
     * @return
     */
    public static Object[] toObjectArray(JSONArray jsonArray) {
        Object[] arr = new Object[jsonArray.length()];
        for (int i = 0; i < arr.length; i ++) {
            arr[i] = jsonArray.get(i);
            if (arr[i] == JSONObject.NULL) // JSON library inserts a non-null marker for NULL values; we need a proper null
                arr[i] = null;
            else {
                try { // see if can be converted to a String
                    arr[i] = ((String) arr[i]).toCharArray()[0];
                } catch (ClassCastException e) {
                }
            }
        }
        return arr;
    }

    public static EnumSeq fromJSON(JSONObject json, Enumerable domain) {
        try {
            String name = json.getString("Name");
            EnumSeq eseq = new EnumSeq.Gappy<>(domain);
            Object[] arr = toObjectArray(json.getJSONArray("Seq"));
            eseq.set(arr);
            eseq.setName(name);
            return eseq;
        } catch (JSONException e) {
            throw new RuntimeException("Error in JSON encoding sequence: " + json.toString());
        }
    }

    public static EnumSeq fromJSON(JSONObject json) {
        try {
            JSONObject jdom = json.getJSONObject("Datatype");
            Enumerable domain = Enumerable.fromJSON(jdom);
            return fromJSON(json, domain);
        } catch (JSONException e) {
            throw new RuntimeException("Error in JSON encoding sequence: " + json.toString());
        }
    }

    public Object[] getStripped() {
        Object[] syms = get();
        int empty = 0;
        for (int i = 0; i < syms.length; i ++) {
            if (syms[i] == null)
                empty += 1;
        }
        Object[] stripped = new Object[syms.length - empty];
        int j = 0;
        for (int i = 0; i < syms.length; i ++) {
            if (syms[i] != null)
                stripped[j ++] = syms[i];
        }
        return stripped;
    }

    public static <T extends dat.Enumerable> List<EnumSeq<T>> loadFasta(String filename, T elementType) throws IOException {
        List<EnumSeq<T>> seqs = new ArrayList<>();
        FastaReader r = new FastaReader(filename, elementType);
        EnumSeq[] rseqs = r.load();
        for (EnumSeq rseq : rseqs) {
            try {
                seqs.add((EnumSeq<T>) rseq);
            } catch (ClassCastException e) {
            }
        }
        r.close();
        return seqs;
    }

    public static EnumSeq parseProtein(String sequence) {
        EnumSeq seq = new EnumSeq(Enumerable.aacid);
        Object[] oarr = new Object[sequence.length()];
        char[] carr = sequence.toCharArray();
        for (int i = 0; i < carr.length; i ++) {
            oarr[i] = (Object) carr[i];
            if (!Enumerable.aacid.isValid(oarr[i]))
                return null;
        }
        seq.set(oarr);
        return seq;
    }

    public static EnumSeq parseDNA(String sequence) {
        EnumSeq seq = new EnumSeq(Enumerable.nacid);
        Object[] oarr = new Object[sequence.length()];
        char[] carr = sequence.toCharArray();
        for (int i = 0; i < carr.length; i ++) {
            oarr[i] = (Object) carr[i];
            if (!Enumerable.nacid.isValid(oarr[i]))
                return null;
        }
        seq.set(oarr);
        return seq;
    }

    public static EnumSeq parseRNA(String sequence) {
        EnumSeq seq = new EnumSeq(Enumerable.nacidRNA);
        Object[] oarr = new Object[sequence.length()];
        char[] carr = sequence.toCharArray();
        for (int i = 0; i < carr.length; i ++) {
            oarr[i] = (Object) carr[i];
            if (!Enumerable.nacidRNA.isValid(oarr[i]))
                return null;
        }
        seq.set(oarr);
        return seq;
    }

    public static EnumSeq<dat.Enumerable> aacid_seq = new EnumSeq<>(dat.Enumerable.aacid);
    public static EnumSeq<dat.Enumerable> nacid_seq = new EnumSeq<>(dat.Enumerable.nacid);

    public static class Gappy<E extends dat.Enumerable> extends EnumSeq<E> {

        public Gappy(E elementType) {
            super(elementType);
        }
        
        @Override
        public boolean isValid(Object value) {
            try {
                Object[] myarr = (Object[]) value;
                for (Object ch : myarr) {
                    if (ch == null)
                        continue;
                    if (!elementType.isValid(ch)) 
                        return false;
                }
            } catch (ClassCastException e) {
                return false;
            }
            return true;
        }

        public static EnumSeq.Gappy fromJSON(JSONObject json, Enumerable domain) {
            try {
                String name = json.getString("Name");
                EnumSeq.Gappy eseq = new EnumSeq.Gappy<>(domain);
                Object[] arr = toObjectArray(json.getJSONArray("Seq"));
                if (!eseq.set(arr))
                    throw new RuntimeException("Invalid sequence in JSON: " + json);
                eseq.setName(name);
                return eseq;
            } catch (JSONException e) {
                throw new RuntimeException("Formatting error in JSON: " + json);
            }
        }

        public static EnumSeq.Gappy fromJSON(JSONObject json) {
            try {
                JSONObject jdom = json.getJSONObject("Datatype");
                Enumerable domain = Enumerable.fromJSON(jdom);
                return fromJSON(json, domain);
            } catch (JSONException e) {
                throw new RuntimeException("Error in JSON encoding sequence: " + json);
            }
        }

        /**
         * New version of the below, takes a buffered reader in, this was to avoid saving files on
         * the GRASP server unnecessarily.
         *
         * @param br
         * @param elementType
         * @param <T>
         * @return
         * @throws IOException
         */
        public static <T extends dat.Enumerable> List<EnumSeq.Gappy<T>> loadClustal(BufferedReader br, T elementType) throws IOException {
            List<EnumSeq.Gappy<T>> seqs = new ArrayList<>();
            AlnReader r = new AlnReader(br, elementType);
            EnumSeq.Gappy[] rseqs = r.load();
            for (EnumSeq.Gappy<T> rseq : rseqs) {
                try {
                    seqs.add((EnumSeq.Gappy<T>) rseq);
                } catch (ClassCastException e) {
                }
            }
            r.close();
            return seqs;
        }

        public static <T extends dat.Enumerable> List<EnumSeq.Gappy<T>> loadClustal(String filename, T elementType) throws IOException {
            List<EnumSeq.Gappy<T>> seqs = new ArrayList<>();
            AlnReader r = new AlnReader(filename, elementType);
            EnumSeq.Gappy[] rseqs = r.load();
            for (EnumSeq.Gappy<T> rseq : rseqs) {
                try {
                    seqs.add((EnumSeq.Gappy<T>) rseq);
                } catch (ClassCastException e) {
                }
            }
            r.close();
            return seqs;
        }

        /**
         * Apparently untested?
         *
         * This is the one used in GRASP atm (08042019) i (ariane) have just repourposed it to use
         * for the buffered reader.
         *
         * @param br buffered reader
         * @param elementType
         * @param gap
         * @param <T>
         * @return
         * @throws IOException
         */
        public static <T extends dat.Enumerable> List<EnumSeq.Gappy<T>> loadFasta(BufferedReader br, T elementType, Character gap) throws IOException {
            //FIXME - untested
            List<EnumSeq.Gappy<T>> seqs = new ArrayList<>();
            FastaReader r = new FastaReader(br, elementType, gap);
            EnumSeq.Gappy[] rseqs = r.loadGappy();
            for (EnumSeq.Gappy rseq : rseqs) {
                try {
                    seqs.add((EnumSeq.Gappy<T>) rseq);
                } catch (ClassCastException e) {
                }
            }
            r.close();
            return seqs;
        }

        public static <T extends dat.Enumerable> List<EnumSeq.Gappy<T>> loadFasta(String filename, T elementType, Character gap) throws IOException {
            //FIXME - untested
            List<EnumSeq.Gappy<T>> seqs = new ArrayList<>();
            FastaReader r = new FastaReader(filename, elementType, gap);
            EnumSeq.Gappy[] rseqs = r.loadGappy();
            for (EnumSeq.Gappy rseq : rseqs) {
                try {
                    seqs.add((EnumSeq.Gappy<T>) rseq);
                } catch (ClassCastException e) {
                }
            }
            r.close();
            return seqs;
        }


    }

    public static Gappy<dat.Enumerable> aacid_gapseq = new Gappy<>(dat.Enumerable.aacid);
    public static Gappy<dat.Enumerable> nacid_gapseq = new Gappy<>(dat.Enumerable.nacid);

    public static class Alignment<E extends dat.Enumerable> {
        
        private final List<EnumSeq.Gappy<E>> seqs;
        private final int width;
        private E domain = null;

        /**
         * Create an alignment structure out of aligned, gappy sequences.
         * @param aseqs 
         */
        public Alignment(List<EnumSeq.Gappy<E>> aseqs) {
            this.seqs = aseqs;
            int w = -1;
            for (EnumSeq.Gappy<E> seq : aseqs) {
                if (w < 0) {
                    w = seq.length();
                    domain = seq.elementType;
                } else if (w != seq.length())
                    throw new RuntimeException("Invalid alignment with sequences of different lengths.");
                else if (domain != seq.elementType)
                    throw new RuntimeException("Invalid alignment with sequences of different types.");
            }
            width = w;
        }

        /**
         * Create an alignment structure out of aligned, gappy sequences.
         * @param aseqs
         */
        public Alignment(EnumSeq.Gappy[] aseqs) {
            this(List.of(aseqs));
        }

        @Override
        public boolean equals(Object other) {
            try {
                Alignment oaln = (Alignment) other;
                int idx = 0;
                for (EnumSeq seq : seqs) {
                    EnumSeq oseq = oaln.getEnumSeq(idx);
                    if (!oseq.equals(seq))
                        return false;
                    if (!oseq.getType().equals(seq.getType()))
                        return false;
                    if (!oseq.getName().equals(seq.getName()))
                        return false;
                    idx += 1;
                }
            } catch (ClassCastException e) {
                return false;
            }
            return true;
        }

        public EnumSeq[] getArray() {
            EnumSeq[] arr = new EnumSeq[seqs.size()];
            seqs.toArray(arr);
            return arr;
        }

        public JSONObject toJSON() {
            JSONObject json = new JSONObject();
            json.put("Datatype", domain.toJSON());
            JSONArray array = new JSONArray();
            for (EnumSeq seq : seqs) {
                JSONObject jseq = seq.toJSONslim(); // no datatype included
                array.put(jseq);
            }
            json.put("Sequences", array);
            return json;
        }

        public static Alignment fromJSON(JSONObject json) {
            try {
                Enumerable domain = Enumerable.fromJSON(json.getJSONObject("Datatype"));
                JSONArray seqarr = json.getJSONArray("Sequences");
                List<EnumSeq> eseqs = new ArrayList<>();
                for (int i = 0; i < seqarr.length(); i++) {
                    JSONObject jseq = seqarr.getJSONObject(i);
                    eseqs.add(EnumSeq.Gappy.fromJSON(jseq, domain));
                }
                Alignment aln = new Alignment(eseqs);
                return aln;
            } catch (JSONException e) {
                throw new RuntimeException("Invalid JSON encoding of alignment");
            }
        }

        public EnumSeq.Gappy<E> getEnumSeq(int index) {
            return seqs.get(index);
        }

        public E getDomain() {
            return domain;
        }

        /**
         * Get the width of alignment--the number of columns.
         * @return number of columns
         */
        public int getWidth() {
            return this.width;
        }
        
        /**
         * Get the height of alignment--the number of sequences.
         * @return number of sequences
         */
        public int getHeight() {
            return seqs.size();
        }
        
        /**
         * Get the names of all the sequences in the alignment
         * @return an array with names
         */
        public String[] getNames() {
            String[] names = new String[seqs.size()];
            for (int i = 0; i < names.length; i ++) {
                EnumSeq.Gappy<E> seq = seqs.get(i);
                names[i] = seq.getName();
            }
            return names;
        }
        
        /**
         * Get the column of enumerable values for a given column, indexed from 0 up to alignment width - 1.
         * @param col column
         * @return array of values in column, null representing gap
         */
        public Object[] getColumn(int col) {
            if (col >= 0 && col < this.width) {
                Object[] syms = new Object[getHeight()];
                for (int i = 0; i < syms.length; i ++)
                    syms[i] = seqs.get(i).get()[col];
                return syms;
            }
            return null;
        }

        /**
         * Get the column with status of gaps (true or false) for a given column, indexed from 0 up to alignment width - 1.
         * @param col column
         * @return array of true or false indicating presence of gap in column
         */
        public Object[] getGapColumn(int col) {
            if (col >= 0 && col < this.width) {
                Object[] syms = new Object[getHeight()];
                for (int i = 0; i < syms.length; i ++)
                    syms[i] = seqs.get(i).get()[col] == null;
                return syms;
            }
            return null;
        }

        /**
         * Determine the number of sequences that have content in this column
         * @param col column
         * @return count of sequences
         */
        public int getOccupancy(int col) {
            if (col >= 0 && col < this.width) {
                int count = 0;
                for (EnumSeq<E> seq : seqs)
                    count += seq.get(col) == null ? 0 : 1;
                return count;
            }
            return 0;
        }
    }
    
}
