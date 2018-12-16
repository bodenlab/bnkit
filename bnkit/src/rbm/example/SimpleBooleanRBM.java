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
import rbm.alg.CD;

import java.io.IOException;

/**
 * This uses the MNIST digits in a Boolean form (grabbed here: http://cs.nyu.edu/%7Eroweis/data.html).
 * The RBM is trained to self-associate digits, i.e. to encode/decode via a hidden layer.
 * @author mikael
 */
public class SimpleBooleanRBM {
    
    public static void main(String[] args) {

        BooleanRBM rbm = new BooleanRBM(28 * 28,20);
        try {
        //    BooleanRBM rbm = new BooleanRBM("/Users/mikael/simhome/rbm/after.tsv");

            EnumVariable[] vars = rbm.getVisibleVars();
            Object[][] data = DataBuf.load("/Users/mikael/simhome/rbm/btrn0.txt", vars, false);
            CD<BooleanRBM> trn = new CD<>(rbm, 1);
            rbm.save("/Users/mikael/simhome/rbm/before.tsv");
            trn.train(data);
            rbm.save("/Users/mikael/simhome/rbm/halfway.tsv");
            trn.train(data);
            rbm.save("/Users/mikael/simhome/rbm/after.tsv");
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
