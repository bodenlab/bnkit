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
import dat.file.TSVFile;
import dat.EnumVariable;
import rbm.BooleanRBM;

import java.io.IOException;
import java.util.Random;

/**
 * This uses the MNIST digits in a Boolean form (grabbed here: http://cs.nyu.edu/%7Eroweis/data.html).
 * The RBM is trained to self-associate digits, i.e. to encode/decode via a hidden layer.
 * @author mikael
 */
public class SimpleBooleanRBM {
    
    public static void main(String[] args) {

        char digit = '1';

        try {
            BooleanRBM rbm = new BooleanRBM("/Users/mikael/simhome/rbm/after.tsv");
            EnumVariable[] vars = rbm.getVisibleVars();
            Object[][] data = DataBuf.load("/Users/mikael/simhome/rbm/btrn" + digit + ".txt", vars, false);
            Object[][] collect0 = new Object[10][];
            Object[][] collect1 = new Object[10][];
            Object[][] collect2 = new Object[10][];
            Object[][] collect3 = new Object[10][];
            for (int i = 0; i < 10; i ++) {
                int n = i;
                Object[] hid = rbm.encode(data[n]);
                Object[] vis = rbm.decode(hid);
                collect0[i] = vis;
            }
            TSVFile tsv = new TSVFile(collect0);
            tsv.save("/Users/mikael/simhome/rbm/bres" + digit + "_0.txt");

            for (int i = 0; i < 10; i ++) {
                int n = i;
                for (int j = 0; j < data[n].length; j ++)
                    if (j > data[n].length / 2)
                        data[n][j] = null;
                Object[] vis = rbm.encode_decode_clamped(data[n]);
                collect1[i] = vis;
            }
            tsv = new TSVFile(collect1);
            tsv.save("/Users/mikael/simhome/rbm/bres" + digit + "_1.txt");

            for (int i = 0; i < 10; i ++) {
                int n = i;
                Object[] vis = rbm.encode_decode_clamped(data[n], 1);
                collect2[i] = vis;
            }
            tsv = new TSVFile(collect2);
            tsv.save("/Users/mikael/simhome/rbm/bres" + digit + "_2.txt");

            for (int i = 0; i < 10; i ++) {
                int n = i;
                Object[] vis = rbm.encode_decode_clamped(data[n], 2);
                collect3[i] = vis;
            }
            tsv = new TSVFile(collect3);
            tsv.save("/Users/mikael/simhome/rbm/bres" + digit + "_3.txt");
            System.exit(0);
/*
            BooleanRBM rbm = new BooleanRBM(28 * 28,20);
            Object[][] data = DataBuf.load("/Users/mikael/simhome/rbm/btrn_all.txt", rbm.getVisibleVars(), false);
            CD<BooleanRBM> trn = new CD<>(rbm, 1);
            rbm.save("/Users/mikael/simhome/rbm/before.tsv");
            trn.train(data);
            rbm.save("/Users/mikael/simhome/rbm/halfway.tsv");
            trn.train(data);
            rbm.save("/Users/mikael/simhome/rbm/after.tsv");
            */
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
