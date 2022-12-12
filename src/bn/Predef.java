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
package bn;

import bn.ctmc.SubstModel;
import bn.ctmc.SubstNode;
import bn.ctmc.matrix.JTT;
import bn.prob.EnumDistrib;
import bn.node.CPT;
import bn.node.CPTPseudo;
import bn.node.DirDT;
import bn.node.GDT;
import dat.Continuous;
import dat.EnumVariable;
import dat.Variable;
import dat.Enumerable;
import java.util.ArrayList;
import java.util.List;

/**
 * This class defines a number of conventional variable and node types.
 * For a variable and node to be specified in a XML file to be loaded or saved, 
 * their types must be defined and mapped in this class. Furthermore, they need
 * to be recognised by methods that check if they accept parameters etc.
 * @author mikael
 */
public class Predef {

    public Predef() {
    }

    public static boolean isParameterised(String typename) {
        if (typename.equalsIgnoreCase("String") || (typename.equalsIgnoreCase("Number")) || typename.equalsIgnoreCase("Distrib")) {
            return true;
        } else {
            return false;
        }
    }

    public static String[] getVariableTypes() {
        return new String[]{"Boolean", "String", "Number", "Real", "Amino acid", "Amino acid extended", "Gap character",
                "Nucleic acid", "Distrib"};
    }

    public static boolean isEnumerable(String typename) {
        if (typename.equalsIgnoreCase("Real") || typename.equalsIgnoreCase("Distrib")) {
            return false;
        } else {
            return true;
        }
    }

    @SuppressWarnings("rawtypes")
    public static Variable getVariable(String varname, String typename, String params) {
        if (typename.equalsIgnoreCase("Boolean")) {
            return Predef.Boolean(varname);
        }
        if (typename.equalsIgnoreCase("String")) {
            String[] values = params.split(";");
            return Predef.Nominal(values, varname);
        }
        if (typename.equalsIgnoreCase("Number")) {
            try {
                Integer max = Integer.parseInt(params);
                return Predef.Number(max, varname);
            } catch (NumberFormatException e) {
                e.printStackTrace();
                return null;
            }
        }
        if (typename.equalsIgnoreCase("Amino acid")) {
            return Predef.AminoAcid(varname);
        }
        if (typename.equalsIgnoreCase("Amino acid extended")) {
            return Predef.AminoAcidExt(varname);
        }
        if (typename.equalsIgnoreCase("Gap character")) {
            return Predef.GapCharacter(varname);
        }
        if (typename.equalsIgnoreCase("Nucleic acid")) {
            return Predef.NucleicAcid(varname);
        }
        if (typename.equalsIgnoreCase("Real")) {
            return Predef.Real(varname);
        }
        if (typename.equalsIgnoreCase("Distrib")) {
            String[] values = params.split(";");
            if (values.length > 1)
                return Predef.Distrib(values, varname);
            else if (!isParameterised(values[0]) && isEnumerable(values[0])) {
                try {
                    Variable temp = getVariable("temp", values[0], null);
                    return Predef.Distrib((Enumerable)temp.getDomain(), varname);
                } catch (RuntimeException ex) {
                    throw new RuntimeException("Invalid specification of variable: " + varname);
                }
            }
        }
        throw new RuntimeException("Invalid specification of variable: " + varname);
    }

    @SuppressWarnings("rawtypes")
    public static BNode getBNode(Variable var, List<Variable> parents, String type) {
        if (parents == null) {
            parents = new ArrayList<>();
        }
        try {
            if (type.equalsIgnoreCase("CPT")) {
                List<EnumVariable> elist = new ArrayList<>();
                for (Variable v : parents) {
                    elist.add((EnumVariable) v);
                }
                return new CPT((EnumVariable) var, elist);
            } else if (type.equalsIgnoreCase("CPTPseudo")) {
                List<EnumVariable> elist = new ArrayList<>();
                for (Variable v : parents) {
                    elist.add((EnumVariable) v);
                }
                return new CPTPseudo((EnumVariable) var, elist);
            }else if (type.equalsIgnoreCase("GDT")) {
                List<EnumVariable> elist = new ArrayList<>();
                for (Variable v : parents) {
                    elist.add((EnumVariable) v);
                }
                return new GDT((Variable<Continuous>) var, elist);
            } else if (type.equalsIgnoreCase("DirDT")) {
                List<EnumVariable> elist = new ArrayList<>();
                for (Variable v : parents) {
                    elist.add((EnumVariable) v);
                }
                return new DirDT(var, elist);
            } else if (type.equalsIgnoreCase("SubstNode")) {
                SubstModel model = new JTT(); // default, to be overridden in "dump" field
                if (parents.size() > 0) {
                    EnumVariable parent = (EnumVariable)parents.get(0);
                    double time = 1.0; // default, to be overridden in "dump" field
                    return new SubstNode((EnumVariable)var, parent, model, time);
                } else {
                    return new SubstNode((EnumVariable)var, model);
                }
            }
            return null;
        } catch (ClassCastException e) {
            return null;
        }
    }

    @SuppressWarnings("rawtypes")
    public static String getType(Variable var) {
        String type = var.getPredef();
        if (type != null) {
            return type;
        }
        return null;
    }

    public static String getType(BNode node) {
        return node.getType();
    }

    /**
     * Get the object representation of a string, for a specified variable.
     *
     * @param var variable (which needs to be a predefined type (@see bn.Variable#setPredef(String))
     * @param vstr string representation of value (as it appears in say a file)
     * @return the Object instance
     */
    public static Object getObject(Variable var, String vstr) {
        try {
            switch (var.getPredef()) {
                case "Boolean": {
                    if (vstr.charAt(0) == '0')
                        return false;
                    else if (vstr.charAt(0) == '1')
                        return true;
                    else
                        return Boolean.parseBoolean(vstr);
                }
                case "Nucleic acid":
                    {
                        Character ch = vstr.charAt(0);
                        if (var.getDomain().isValid(ch)) {
                            return ch;
                        }       break;
                    }
                case "Amino acid":
                    {
                        Character ch = vstr.charAt(0);
                        if (var.getDomain().isValid(ch)) {
                            return ch;
                        }       break;
                    }
                case "Amino acid extended":
                {
                    Character ch = vstr.charAt(0);
                    if (var.getDomain().isValid(ch)) {
                        return ch;
                    }       break;
                }
                case "Gap character":
                {
                    Character ch = vstr.charAt(0);
                    if (var.getDomain().isValid(ch)) {
                        return ch;
                    }       break;
                }
                case "String":
                    if (var.getDomain().isValid(vstr)) {
                        return vstr;
                    }   break;
                case "Number":
                    Integer n = Integer.parseInt(vstr);
                    if (var.getDomain().isValid(n)) {
                        return n;
                    }   break;
                case "Real":
                    Double y = Double.parseDouble(vstr);
                    if (var.getDomain().isValid(y)) {
                        return y;
                    }   break;
                case "Distrib":
                    EnumDistrib d = EnumDistrib.parseEnumDistrib(vstr, (Enumerable)var.getDomain());
                    if (var.getDomain().isValid(d)) {
                        return d;
                    }   break;
            }
        } catch (NumberFormatException e) {
            return null;
        }
        return null;
    }

    public static EnumVariable Boolean(String name) {
        EnumVariable var = new EnumVariable(Enumerable.bool, name);
        var.setPredef("Boolean");
        return var;
    }

    public static EnumVariable Boolean() {
        return Predef.Boolean("Bool");
    }

    public static EnumVariable NucleicAcid(String name) {
        EnumVariable var = new EnumVariable(Enumerable.nacid, name);
        var.setPredef("Nucleic acid");
        return var;
    }

    public static EnumVariable NucleicAcid() {
        return Predef.NucleicAcid("NT");
    }

    public static EnumVariable AminoAcid(String name) {
        EnumVariable var = new EnumVariable(Enumerable.aacid, name);
        var.setPredef("Amino acid");
        StringBuilder sbuf = new StringBuilder("");
        for (Object v : Enumerable.aacid.getValues()) {
            sbuf.append(v).append(";");
        }
        var.setParams(sbuf.toString());
        return var;
    }
    public static EnumVariable AminoAcidExt(String name) {
        EnumVariable var = new EnumVariable(Enumerable.aacid_ext, name);
        var.setPredef("Amino acid extended");
        StringBuilder sbuf = new StringBuilder("");
        for (Object v : Enumerable.aacid_ext.getValues()) {
            sbuf.append(v).append(";");
        }
        var.setParams(sbuf.toString());
        return var;
    }
    public static EnumVariable GapCharacter(String name) {
        EnumVariable var = new EnumVariable(Enumerable.gap_character, name);
        var.setPredef("Gap character");
        StringBuilder sbuf = new StringBuilder("");
        for (Object v : Enumerable.gap_character.getValues()) {
            sbuf.append(v).append(";");
        }
        var.setParams(sbuf.toString());
        return var;
    }

    public static EnumVariable AminoAcid() {
        return Predef.AminoAcid("AA");
    }
    
    public static EnumVariable SecondaryStructure(String name) {
    	EnumVariable var = new EnumVariable(Enumerable.secondaryStructure, name);
    	var.setPredef("Secondary Structure");
    	var.setParams("H;E;C;");
    	return var;
    }

    public static EnumVariable Nominal(String[] values, String name) {
        Enumerable domain = new Enumerable(values);
        EnumVariable var = new EnumVariable(domain, name);
        StringBuilder sbuf = new StringBuilder("");
        for (Object v : values) {
            sbuf.append(v).append(";");
        }
        var.setParams(sbuf.toString());
        var.setPredef("String");
        return var;
    }

    public static EnumVariable Nominal(String... values) {
        return Predef.Nominal(values, "Nom");
    }

    public static EnumVariable Number(int max, String name) {
        Enumerable domain = new Enumerable(max);
        EnumVariable var = new EnumVariable(domain, name);
//        var.setName("Number");
        var.setParams("" + max);
        var.setPredef("Number");
        return var;
    }

    public static EnumVariable Number(int max) {
        return Predef.Number(max, "Num");
    }

    public static Variable<Continuous> Real(String name) {
        Continuous domain = new Continuous();
        Variable<Continuous> var;
        var = new Variable<>(domain, name);
        var.setPredef("Real");
        return var;
    }

    public static Variable<Continuous> Real() {
        return Predef.Real("Rl");
    }

    public static Variable<EnumDistrib> Distrib(String... values) {
        return Predef.Distrib(values, "Distrib");
    }
    
    public static Variable<EnumDistrib> Distrib(Enumerable dom, String name) {
        EnumDistrib domain = new EnumDistrib(dom);
        Variable<EnumDistrib> var = new Variable<>(domain, name);
        var.setPredef("Distrib");
        StringBuilder sbuf = new StringBuilder("");
        for (Object v : dom.getValues()) {
            sbuf.append(v).append(";");
        }
        var.setParams(sbuf.toString());
        return var;
    } 
    
    public static Variable<EnumDistrib> Distrib(String[] values, String name) {
        return Predef.Distrib(new Enumerable(values), name);
    } 
}
