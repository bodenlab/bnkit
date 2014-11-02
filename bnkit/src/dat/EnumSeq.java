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
    
    public EnumSeq(E elementType) {
        super(elementType);
    }
    
    public void setName(String name) {
        this.name = name;
    }
    
    public String getName() {
        return name;
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
            for (EnumSeq rseq : rseqs) {
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
        
        public Alignment(List<EnumSeq.Gappy<E>> aseqs) {
            
        }
        
    }
    
}
