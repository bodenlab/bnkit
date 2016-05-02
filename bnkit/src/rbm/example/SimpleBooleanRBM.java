/*
    bnkit -- software for building and using probabilistic models
    including Bayesian networks.
    Copyright (C) 2014-2016  M. Boden et al.

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
package rbm.example;

import bn.file.DataBuf;
import dat.EnumVariable;
import rbm.BooleanRBM;

/**
 *
 * @author mikael
 */
public class SimpleBooleanRBM {
    
    static BooleanRBM rbm = new BooleanRBM(5,2);

    public static void main(String[] args) {
        
        EnumVariable[] vars = rbm.getVisibleVars();
        Object[][] data = DataBuf.load("/Users/mikael/simhome/rbm/example.txt", vars, false);
        for (int c = 0; c < data[0].length; c ++) {
            System.out.print(vars[c].getName() + "\t");
        }
        System.out.println();
        for (int r = 0; r < data.length; r ++) {
            for (int c = 0; c < data[r].length; c ++) {
                System.out.print(data[r][c] + "\t");
            }
            System.out.print("==>\t");
            Object[] hid = rbm.encode(data[r]);
            for (int h = 0; h < hid.length; h ++) {
                System.out.print(hid[h] + "\t");
            }
            System.out.print("==>\t");
            Object[] recoded = rbm.decode(hid);
            for (int i = 0; i < recoded.length; i ++) {
                System.out.print(recoded[i] + "\t");
            }
            System.out.println();
            
        }
        
    }
}
