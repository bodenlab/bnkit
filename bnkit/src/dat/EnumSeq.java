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
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author mikael
 * @param <E>
 */
public class EnumSeq<E extends Enumerable> extends SeqDomain<E> {

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
    
    @Override
    public E getType() {
        return super.getType();
    }
    
    public String getInfo() {
    	return info;
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
    
    public static <T extends Enumerable> List<EnumSeq<T>> loadFasta(String filename, T elementType) throws IOException {
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
    
    public static EnumSeq<Enumerable> aacid_seq = new EnumSeq<>(Enumerable.aacid);
    public static EnumSeq<Enumerable> nacid_seq = new EnumSeq<>(Enumerable.nacid);

    public static class Gappy<E extends Enumerable> extends EnumSeq<E> {

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

        public static <T extends Enumerable> List<EnumSeq.Gappy<T>> loadClustal(String filename, T elementType) throws IOException {
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

        public static <T extends Enumerable> List<EnumSeq.Gappy<T>> loadFasta(String filename, T elementType, Character gap) throws IOException {
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

    public static Gappy<Enumerable> aacid_gapseq = new Gappy<>(Enumerable.aacid);
    public static Gappy<Enumerable> nacid_gapseq = new Gappy<>(Enumerable.nacid);

    public static class Alignment<E extends Enumerable> {
        
        private final List<EnumSeq.Gappy<E>> seqs;
        private final int width;
        
        /**
         * Create an alignment structure out of aligned, gappy sequences.
         * @param aseqs 
         */
        public Alignment(List<EnumSeq.Gappy<E>> aseqs) {
            this.seqs = aseqs;
            int w = -1;
            for (EnumSeq.Gappy<E> seq : aseqs) {
                if (w < 0)
                    w = seq.length();
                else if (w != seq.length())
                    throw new RuntimeException("Invalid alignment with sequences of different lengths.");
            }
            width = w;
        }
        
        public EnumSeq.Gappy<E> getEnumSeq(int index) {
            return seqs.get(index);
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
    }
    
}
