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
package bn.file;

import dat.Variable;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import bn.*;

/**
 * Class that defines methods for loading and saving data for training and
 * testing Bayesian networks.
 *
 * @author mikael
 */
public class DataBuf {

    /**
     * Load a data file containing tab-separated values for variables. Return
     * values for listed variables only, and in the specified order. If the
     * variable is not available (according to the header) then its values are
     * set to null. Note that variables need to be of pre-defined types (@see bn.Predef).
     * @param filename
     * @param nodes the variables (as nodes in a Bayesian network)
     * @return all values in their native format [row][field]
     */
    public static Object[][] load(String filename, List<BNode> nodes) {
        return load(filename, nodes, true);
    }

    /**
     * Load a data file containing tab-separated values for variables. Return
     * values for listed variables only, and in the specified order. If the
     * variable is not available (according to an optional header) then its
     * values are set to null. If no header is used, then values are assumed to
     * be ordered as the nodes. Note that variables need to be of pre-defined
     * types (@see bn.Predef).
     * @param filename
     * @param nodes the variables (as nodes in a Bayesian network)
     * @param hasHeader true if the file contains a header listing available
     * variables (by name), false if no header is included
     * @return all values in their native format [row][field]
     */
    public static Object[][] load(String filename, List<BNode> nodes, boolean hasHeader) {
        Variable[] var_arr = new Variable[nodes.size()];
        for (int i = 0; i < nodes.size(); i++) {
            var_arr[i] = nodes.get(i).getVariable();
        }
        return DataBuf.load(filename, var_arr, hasHeader);
    }

    /**
     * Load a data file containing tab-separated values for variables. Return
     * values for listed variables only, and in the specified order. If the
     * variable is not available (according to the header) then its values are
     * set to null. Note that variables need to be of pre-defined types (@see bn.Predef).
     * @param filename
     * @param vars the variables
     * @return all values in their native format [row][field]
     */
    public static Object[][] load(String filename, Variable[] vars) {
        return load(filename, vars, true);
    }

    /**
     * Load a data file containing tab-separated values for variables. Return
     * values for listed variables only, and in the specified order. If the
     * variable is not available (according to an optional header) then its
     * values are set to null. If no header is used, then values are assumed to
     * be ordered as the nodes. Note that variables need to be of pre-defined
     * types (@see bn.Predef).
     * @param filename
     * @param vars the variables
     * @param hasHeader true if the file contains a header listing available
     * variables (by name), false if no header is included
     * @return all values in their native format [row][field]
     */
    public static Object[][] load(String filename, Variable[] vars, boolean hasHeader) {
        File file = new File(filename);
        int nvars_list = vars.length;
        int nvars_file = -1;
        try {
            BufferedReader reader = new BufferedReader(new FileReader(file));
            int[] field_index = new int[vars.length];
            List<Object[]> lines = new ArrayList<Object[]>();
            String line = null;
            try {
                int cnt = 1;
                while ((line = reader.readLine()) != null) {
                    String[] words = line.trim().split("\t");
                    if (nvars_file < 0) {
                        nvars_file = words.length;
                    } else if (nvars_file != words.length) {
                        throw new RuntimeException("Invalid number of values in \"" + line + "\"");
                    }
                    if (cnt == 1 && hasHeader) { // first row is header
                        for (int j = 0; j < nvars_list; j++) {
                            boolean found = false; // has not found a matching header for variable
                            for (int i = 0; i < words.length; i++) {
                                if (words[i].equals(vars[j].getName())) {
                                    field_index[j] = i;
                                    found = true;
                                    break;
                                }
                            }
                            if (!found) // no variable for this header, ignore later
                            {
                                field_index[j] = -1;
                            }
                        }
                    } else {
                        if (cnt == 1) { // no header, prepare fictitious header map
                            for (int j = 0; j < nvars_list; j++) {
                                field_index[j] = j;
                            }
                        }
                        Object[] values = new Object[nvars_list];
                        for (int i = 0; i < nvars_list; i++) {
                            Variable var = vars[i];
                            if (field_index[i] == -1) {
                                values[i] = null;
                            } else {
                                String vstr = words[field_index[i]];
                                if (vstr.equalsIgnoreCase("null") || vstr.equalsIgnoreCase("nil") || vstr.equals("-")) {
                                    values[i] = null;
                                } else {
                                    values[i] = Predef.getObject(var, vstr);
                                    if (values[i] == null) {
                                        throw new RuntimeException("Invalid value in \"" + line + "\" (field " + (i + 1) + " for variable " + var.getName() + " of type " + var.getPredef() + ")");
                                    }
                                }
                            }
                        }
                        lines.add(values);
                    }
                    cnt++;
                }
                Object[][] all = new Object[lines.size()][vars.length];
                for (int row = 0; row < lines.size(); row++) {
                    Object[] values = lines.get(row);
                    for (int col = 0; col < values.length; col++) {
                        all[row][col] = values[col];
                    }
                }
                return all;
            } catch (IOException e) {
                e.printStackTrace();
                return null;
            }
        } catch (FileNotFoundException e1) {
            e1.printStackTrace();
            return null;
        }
    }

    public static void main(String[] argv) {
        // Group	Alpha	Cpl	Age	Gender
        Variable v1 = Predef.Nominal(new String[]{"healthyBP", "healthyMDD", "atrisk", "bipolar", "mel", "nonmel"}, "Class");
        Variable v2 = Predef.Real("Alpha");
        Variable v3 = Predef.Real("Cpl");
        Variable v4 = Predef.Real("Age");
        Variable v5 = Predef.Nominal(new String[]{"Male", "Female"}, "Gender");
        Object[][] data = DataBuf.load("data/short.cas", new Variable[]{v1, v2, v3, v5, v4});
        for (Object[] row : data) {
            for (int f = 0; f < row.length; f++) {
                System.out.print(row[f] + "\t");
            }
            System.out.println();
        }
    }
}
