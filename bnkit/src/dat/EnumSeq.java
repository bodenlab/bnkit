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

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author mikael
 * @param <E>
 */
public class EnumSeq<E extends Enumerable> extends SeqDomain<E> {

    public EnumSeq(E elementType) {
        super(elementType);
    }
    
    public static <T extends Enumerable> List<EnumSeq<T>> loadFasta(String filename) {
        List<EnumSeq<T>> seqs = new ArrayList<>();
        // ... TODO
        
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
                for (Object myarr1 : myarr) {
                    if (myarr1 == null)
                        continue;
                    if (!elementType.isValid(myarr1)) 
                        return false;
                }
            } catch (ClassCastException e) {
                return false;
            }
            return true;
        }

        public static <T extends Enumerable> List<EnumSeq.Gappy<T>> loadClustal(String filename) {
            List<EnumSeq.Gappy<T>> seqs = new ArrayList<>();
            // ... TODO

            return seqs;
        }
    }

    public static Gappy<Enumerable> aacid_gapseq = new Gappy<>(Enumerable.aacid);
    public static Gappy<Enumerable> nacid_gapseq = new Gappy<>(Enumerable.nacid);
}
